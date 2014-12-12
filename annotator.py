#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import json
import mwparserfromhell

from collections import defaultdict, OrderedDict
from cruzdb import Genome
from wikitools import wiki, category, page


def parse_args():
    """
    Returns the arguments from the command line
    """

    # create argparse object with description of script
    program_description = "Phenotypically profile an individual using variant information."
    parser = argparse.ArgumentParser(description=program_description)

    # positional arguments
    parser.add_argument("-f", type=str, required=True,
                        dest="vcf_file", help="File containing variants in .vcf format")
    parser.add_argument("-v", dest="verbose", action="store_true", help="Verbose mode")

    # parse arguments and return list
    return parser.parse_args()


class variant():
    """
    Holds all info regarding variant
    """

    def __init__(self, param):
        self.CHROM = "chr" + param[0]
        self.POS = int(param[1])
        self.ID = param[2]
        self.REF = param[3]
        self.ALT = param[4].split(",")
        #self.QUAL = param[5]
        self.FILTER = param[6]
        self.GENOTYPE = "(" + self.get_genotype(param[8]) + ")"


    def __str__(self):
        return self.ID + self.GENOTYPE

    def get_genotype(self, data):
        """
        Returns the genotype based on the format in the vcf
        """
        genotype = data.split(":")[0].split("/")
        if genotype[0] == ".": return "unknown"     # check for undetermined genotype

        # build dict to map code from vcf to allelle
        code2allele = {str(index+1): a for index, a in enumerate(self.ALT)}
        code2allele["0"] = self.REF

        return ";".join([code2allele[gt] for gt in genotype])


def gene2omim(morbidmap):
    """
    Returns a dict mapping gene/locus symbols to MIM numbers
    """

    genemap = defaultdict(dict)

    # parse morbidmap
    with open(morbidmap, "rb") as f:
        lines = [line.strip().split('|') for line in f.readlines()]

    ignore_chars = ["?", "[", "{"]      # bad identifications

    for line in lines:
        if line[0][0] in ignore_chars:
            continue
        for gene in line[1].split(','):
            # grab diagnosis, omim number and cytogenetic location
            genemap[gene.strip()]["omim"] = line[2]
            genemap[gene.strip()]["dx"] = line[0]
            genemap[gene.strip()]["loc"] = line[3]

    return genemap


def omim_parser(omim):
    """
    Returns a dict mapping OMIM number to variant information
    """

    omimmap = defaultdict(dict)

    # fields we care about capturing from the OMIM file
    fieldcapture = ["NO", "TX", "CS"]

    # parse omim file into string
    with open(omim, "rb") as f:
        omim = "\n".join(line.strip() for line in f.readlines())

    # split the string to each record in the OMIM file
    records = omim.split("*RECORD*")

    for index, record in enumerate(records):
        # zeroth index is blank, skip it
        if index == 0:
            continue

        field_dict = dict()
        fields = record.split("*FIELD* ")

        for param in fields:
            subfields = param.split("\n")
            # head will be NO, TI, TX, AV, etc.
            head = subfields[0]

            # ignore fields we don't care about
            if head not in fieldcapture:
                continue

            # grab all of the text in the field
            summary = "\n".join(subfields[1:]).strip()

            # only want description subfield in summary
            if head == "TX":
                text = summary.split("\n\n")
                try:
                    desc_index = text.index("DESCRIPTION") + 1
                except:
                    # record does not contain DESCRIPTION
                    continue

                field_dict[head] = text[desc_index]
            else:
                field_dict[head] = summary

        try:
            # some entries do not have a "CS" section do a try to catch KeyError
            omimmap[field_dict["NO"]]["TX"] = field_dict["TX"]
            omimmap[field_dict["NO"]]["CS"] = field_dict["CS"]
        except(KeyError):
            pass

    return omimmap


def pos2gene(chr, pos, g):
    """
    Returns a gene name for a given chromosome and position
    chr type(str): "chr2"
    pos type(int): 14325 
    """

    # find the nearest neighbor to the passed in position
    try:
        nearest = g.knearest("refGene", chr, pos, pos+1, k=1)
    except:
        return "non-coding"

    # grab first (closest) gene
    gene = nearest[0]
    exon_list = gene.exons

    # check for position in exon regions return immediately if found
    for beg, end in exon_list:
        if pos >= beg and pos <= end:
            return gene.name2
    # if here the position is non-coding
    return "non-coding"


def mutation_search(variant, genome, unique):
    """
    For variants that do not have an rsID number check if position
    falls in the coding region of a gene
    """

    # get gene symbol and add it to the set if it is coding
    gene = pos2gene(variant.CHROM, variant.POS, g)
    if gene != "non-coding":
        unique_genes[gene] += 1


def snp_search(variant, site, snp_output):
    """
    For variants with an rsID look up the rsID + Genotype in SNPedia
    for diagnostic information
    """

    data = {}

    # fields we want to capture in snpedia
    fields = set(["magnitude", "repute", "summary"])

    # grab the page for the snp
    snp = variant.ID + variant.GENOTYPE
    pagehandle = page.Page(site,snp)

    # check if page exists, else return
    try:
        snp_page = pagehandle.getWikiText()
    except:
        # page not found for snp
        return

    # grab all of the data from the page
    wikicode = mwparserfromhell.parse(snp_page)
    template = wikicode.filter_templates()[0]

    if template.name.strip() == "Genotype":     # check template is named correctly
        for param in template.params:
            p = param.strip().split("=")
            if p[0] in fields:      # grab only the params we care about
                data[p[0]] = p[1]

    if not data:
        return
    else:
        pos = str(variant.CHROM) + ":" + str(variant.POS)
        try:
            # ignore common snps
            if data["summary"].lower().startswith("common"):
                return
            snp_output[variant.ID].append(("Position", pos))
            snp_output[variant.ID].append(("Genotype",variant.GENOTYPE))
            snp_output[variant.ID].append(("Phenotype", data["summary"]))
            snp_output[variant.ID].append(("Repute", data["repute"]))
            snp_output[variant.ID].append(("Importance", data["magnitude"]))
        except:
            # some entires are incomplete remove and continue
            snp_output.pop(variant.ID, None)


def vcf_parser(vcf_file, g, unique_genes, snp_output, site, sex, v):
    """
    Parse the passed in vcf file
    """

    with open(vcf_file, "rb") as vcf:
        for index, line in enumerate(vcf):
            if line.startswith("#"):    # ignore header info
                continue
            if v:
                if index % 10000 == 0:
                    print "%d variants complete" % index

            # read the line into a variant object
            v = variant(line.strip().split("\t"))

            # filter out low quality variants
            if v.FILTER != "PASS" or v.GENOTYPE == "(unknown)":
                continue

            if v.CHROM.lower() == "y":
                sex = "Male"

            # for non-SNPs find the gene and add it to the list
            if v.ID == ".":
                mutation_search(v, g, unique_genes)
            elif v.ID in known_snps:
                snp_search(v, site, snp_output)
            else:
                continue



def gene_mapping(genes, genemap, mut_output):
    """
    Build mutation object using gene symbols and OMIM queries
    """

    for gene, num in genes.iteritems():

        try:
            # get OMIM information
            omimid = genemap[gene]["omim"]
            dx = genemap[gene]["dx"]
            loc = genemap[gene]["loc"]
        except(KeyError):
            # does not correspond to a known OMIM ID
            continue

        entry = omimmap[omimid]

        if len(entry) == 0:
            continue

        try:
            # build object
            mut_output[gene].append(("Diagnosis",dx))
            mut_output[gene].append(("Description",entry["TX"]))
            mut_output[gene].append(("Cytogenetic location",loc))
            mut_output[gene].append(("Number of Mutations",num))
            mut_output[gene].append(("Synoposis",entry["CS"]))
        except(KeyError):   # sometimes no Synopsis. Catch KeyError
            continue


def write_data(mutations, snpedia):
    """
    Write json objects to disk to be read by poseidon
    """

    print "\nWriting 'mutations.json', and 'snpedia.json'. Upload to Poseidon for viewing."
    with open("mutations.json", "wb") as outfile:
        json.dump(mutations, outfile, indent=4)
    with open("snpedia.json", "wb") as outfile:
        json.dump(snpedia, outfile, indent=4)

if __name__ == '__main__':

    # arguments from command line
    args = parse_args()

    # grab command line arguments
    vcf_file = args.vcf_file
    v = args.verbose

    if v:
        print "\nLoading local genome and list of snps in snpedia...\n"
    try:
        g = Genome("sqlite:///db/hg19.db")
        with open("db/known_snps.txt", "rb") as f:
            known_snps = set([snp.strip() for snp in f])
    except:
        print "Local data not found. Please run pre_annotator.py first."
        exit(1)

    if v:
        print "\nMapping OMIM data...\n"
    # get a dict mapping gene symbols to OMIM IDs
    genemap = gene2omim("omim/morbidmap")
    # get a dict mapping OMIM IDs to OMIM objects
    omimmap = omim_parser("omim/omim.txt")

    # site for querying snps
    site = wiki.Wiki("http://bots.snpedia.com/api.php")

    # create object to store output for mutations and for snps
    mut_output = defaultdict(list)
    snp_output = defaultdict(list)

    # sets to hold unique genes for quick parsing of same data after an initial run
    #unique_genes = set()
    unique_genes = defaultdict(int)

    # initialize sex to be female
    sex = "Female"

    if v:
        print "\nBegin parsing vcf by querying SNPedia and finding gene mutations...\n"
    vcf_parser(vcf_file, g, unique_genes, snp_output, site, sex, v)


    if v:
        print "\nMapping gene symbols to OMIM entries...\n"
    gene_mapping(unique_genes, genemap, mut_output)

    # sort on number of mutations
    mut_sorted = OrderedDict(sorted(mut_output.iteritems(), key=lambda (k,v): v[3][1],reverse=True))

    # insert individual descriptor keys to each object
    mut_sorted["Sex"] = sex
    mut_sorted["Desc"] = ("The following is a list of genes and a possible diagnosis if the gene were "
                          "knocked out. For all mutations without an associated SNP, there is a check "
                          "to see if the mutation position is in the coding region of a gene. The number "
                          "of mutations in each gene is given and the genes are sorted in decreasing "
                          "number of mutations. It is unknown if the given mutations are deleterious.")
    snp_output["Sex"] = sex
    snp_output["Desc"] = ("The following is a list of SNPs found in SNPedia and their associated "
                          "phenotypes given the invidual's genotype. Repute determines whether the "
                          "phenotype is 'Good' or 'Bad'. Importance is a community derived number "
                          "of how interesting/important the SNP is (scale: 0-10). Common SNPs are "
                          "not shown.")

    # write json objects
    write_data(mut_sorted, snp_output)
