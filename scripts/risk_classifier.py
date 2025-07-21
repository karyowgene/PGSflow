if genotype == "CC":
    risk = OR ** 2  # Homozygous risk allele
elif genotype == "CT":
    risk = OR ** 1  # Heterozygous
else:
    risk = 1.0      # No risk allele
