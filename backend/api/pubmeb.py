import requests
import xml.etree.ElementTree as ET

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PUBMED_URL = "https://pubmed.ncbi.nlm.nih.gov/"

# Search PubMed for "Dental AI"
params = {
    "db": "pubmed",
    "term": "Management of Crohn's Disease in the Dental Setting",
    "retmode": "json",
    "retmax": 5  # Number of results
}

response = requests.get(BASE_URL + "esearch.fcgi", params=params)
data = response.json()

# Get the list of PubMed IDs (PMIDs)
pmids = data["esearchresult"]["idlist"]
print("PubMed IDs:", pmids)

def fetch_article_details(pmids):
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),  # Pass comma-separated PMIDs
        "retmode": "xml"  # Use XML instead of JSON
    }
    
    response = requests.get(BASE_URL + "efetch.fcgi", params=params)
    
    # Parse XML response
    root = ET.fromstring(response.text)
    
    articles = []
    for article, pmid in zip(root.findall(".//PubmedArticle"), pmids):
        title = article.find(".//ArticleTitle").text if article.find(".//ArticleTitle") is not None else "No Title"
        abstract = article.find(".//AbstractText").text if article.find(".//AbstractText") is not None else "No Abstract"
        article_link = f"{PUBMED_URL}{pmid}/"  # Construct the article link
        articles.append({"title": title, "abstract": abstract, "link": article_link})
    
    return articles

# Fetch article details
article_data = fetch_article_details(pmids)
print(article_data)
# Print results
for i, article in enumerate(article_data, 1):
    print(f"\nArticle {i}:")
    print(f"Title: {article['title']}")
    print(f"Abstract: {article['abstract']}...")  # Print only the first 300 characters
    print(f"Link: {article['link']}")
