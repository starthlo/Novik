from Bio import Entrez
from difflib import SequenceMatcher

Entrez.email = "chatgptsuperpowers@gmail.com"

def similarity(a, b):
    return SequenceMatcher(None, a.lower(), b.lower()).ratio()

def get_pubmed_links(title, max_results=3):
    # Search PubMed
    handle = Entrez.esearch(db="pubmed", term=title, retmax=10)
    record = Entrez.read(handle)
    ids = record["IdList"]

    if not ids:
        return []

    # Fetch articles
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.read(handle)

    # Collect titles and similarity
    matches = []
    for article in records["PubmedArticle"]:
        article_title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        pmid = article["MedlineCitation"]["PMID"]
        sim_score = similarity(title, article_title)
        matches.append((sim_score, article_title, f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"))

    # Return top N by similarity
    top_matches = sorted(matches, key=lambda x: x[0], reverse=True)[:max_results]
    return top_matches

# Example titles
titles = [
    "Managing Patients with Dental Procedures while on Immunosuppressants",
    "Management of Crohn's Disease in the Dental Setting",
    "Immunosuppressive Therapy and Outcomes"
]

# Print results
for t in titles:
    print(f"\nTop matches for: {t}")
    links = get_pubmed_links(t)
    for i, (score, found_title, url) in enumerate(links, 1):
        print(f"{i}. {found_title} - {url}")
