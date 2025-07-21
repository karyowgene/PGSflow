import requests
import pandas as pd

# Fetch BMD-associated SNPs from GWAS Catalog API
def fetch_gwas_snps(trait, pval):
    url = "https://www.ebi.ac.uk/gwas/rest/api/associations"
    params = {
        "efoTrait": trait,
        "pValue": f"<{pval}",
        "size": 1000
    }
    response = requests.get(url, params=params)
    data = response.json()["_embedded"]["associations"]
    
    # Extract SNP info
    snps = []
    for assoc in data:
        for snp in assoc["loci"]:
            snps.append({
                "rsid": snp["strongestRiskAlleles"][0]["riskAlleleName"].split("-")[0],
                "effect_allele": snp["strongestRiskAlleles"][0]["riskAllele"],
                "or": snp["orValue"],
                "pvalue": assoc["pvalue"]
            })
    return pd.DataFrame(snps)

if __name__ == "__main__":
    df = fetch_gwas_snps(snakemake.params["trait"], snakemake.params["pval"])
    df.to_csv(snakemake.output["snps"], sep="\t", index=False)
