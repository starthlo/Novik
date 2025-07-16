import os
from typing import Dict, List

from Bio import Entrez
from django.utils import timezone
from dotenv import load_dotenv

load_dotenv()

Entrez.email = os.environ.get("ENTREZ_EMAIL")


def get_articles(term: str, count: int = 3) -> List[Dict[str, str]]:
    ymax = timezone.now().date().year
    ymin = ymax - 5

    handle = Entrez.esearch(
        db="pubmed",
        term=term,
        datetype="pdat",
        mindate=ymin,
        maxdate=ymax,
        retmax=count,
        sort="relevance",
    )
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return []

    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    articles = []
    for article in records["PubmedArticle"]:
        art = article["MedlineCitation"]["Article"]
        title = art["ArticleTitle"]
        pmid = article["MedlineCitation"]["PMID"]

        authors = art.get("AuthorList", [])
        author_names = []
        for a in authors[:6]:
            last = a.get("LastName", "")
            initials = a.get("Initials", "")
            if last:
                author_names.append(f"{last} {initials}")

        if len(authors) > 6:
            author_names.append("et al.")

        author_str = ", ".join(author_names)

        journal = art["Journal"]["Title"]
        year = art["Journal"]["JournalIssue"]["PubDate"].get("Year", "n.d.")
        volume = art["Journal"]["JournalIssue"].get("Volume", "")
        issue = art["Journal"]["JournalIssue"].get("Issue", "")
        pages = art.get("Pagination", {}).get("MedlinePgn", "")

        dois = [
            str(el)
            for el in art.get("ELocationID", [])
            if getattr(el, "attributes", {}).get("EIdType") == "doi"
        ]
        doi = f"doi:{dois[0]}" if dois else ""
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

        articles.append(
            {
                "title": title,
                "authors": author_str,
                "journal": journal,
                "year": year,
                "volume": volume,
                "issue": issue,
                "pages": pages,
                "doi": doi,
                "url": url,
            }
        )

    return articles
