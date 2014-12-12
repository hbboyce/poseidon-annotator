#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import mwparserfromhell
import os

from cruzdb import Genome
from wikitools import wiki, category, page

directory = "db"

print "Creating db/ directory"
if not os.path.exists(directory):
    os.makedirs(directory)

print "Saving local version of hg19"
Genome('hg19').mirror(["refGene"], "sqlite:///db/hg19.db")

print "Saving list of SNPs in SNPedia"
site = wiki.Wiki("http://bots.snpedia.com/api.php")
known_snps = set()
snps = category.Category(site, "Is_a_snp")
for article in snps.getAllMembersGen(namespaces=[0]):
    if article.title.lower().startswith("rs"):
        known_snps.add(article.title.lower())

with open("db/known_snps.txt", "wb") as f:
    for snp in known_snps:
        f.write(snp + "\n")
