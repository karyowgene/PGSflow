import requests

# Fetch PRS weights from PGS Catalog
def fetch_pgs_weights(pgs_id):
    url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
    response = requests.get(url)
    data = response.json()
    return data["variants"]

if __name__ == "__main__":
    weights = fetch_pgs_weights(snakemake.params["pgs_id"])
    with open(snakemake.output["weights"], "w") as f:
        f.write("rsid\tweight\teffect_allele\n")
        for variant in weights:
            f.write(f"{variant['rsID']}\t{variant['weight']}\t{variant['effect_allele']}\n")
