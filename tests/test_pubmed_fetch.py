import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from src.literature.pubmed_fetch import fetch_pubmed_data

@pytest.fixture
def mock_config():
    return {
        "entrez_email": "test@example.com",
        "journals": ["Nature", "Science"],
        "keywords": ["CRISPR", "Gene"],
        "max_records": 10
    }

@pytest.fixture
def sample_esearch_result():
    """Symuluje odpowiedź z listą ID artykułów."""
    return {"IdList": ["12345", "67890"]}

@pytest.fixture
def sample_efetch_result():
    """
    Symuluje odpowiedź XML z detalami artykułów (po sparsowaniu przez Entrez.read).
    """
    # Pełne dane z DOI
    article1 = {
        "MedlineCitation": {
            "PMID": "12345",
            "Article": {
                "ArticleTitle": "CRISPR advances",
                "Journal": {"Title": "Nature"}
            }
        },
        "PubmedData": {
            "ArticleIdList": [
                MagicMock(attributes={"IdType": "pubmed"}, __str__=lambda x: "12345"),
                MagicMock(attributes={"IdType": "doi"}, __str__=lambda x: "10.1038/nature123")
            ]
        }
    }

    # Brak DOI
    article2 = {
        "MedlineCitation": {
            "PMID": "67890",
            "Article": {
                "ArticleTitle": "Gene editing basics",
                "Journal": {"Title": "Science"}
            }
        },
        "PubmedData": {
            "ArticleIdList": [
                MagicMock(attributes={"IdType": "pubmed"}, __str__=lambda x: "67890")
            ]
        }
    }
    
    return {"PubmedArticle": [article1, article2]}

# ---- Testy ----

@patch("Bio.Entrez.efetch")
@patch("Bio.Entrez.esearch")
def test_parsing_logic(mock_esearch, mock_efetch, mock_config, sample_esearch_result, sample_efetch_result):

    # Mocki danych
    mock_esearch.return_value.__enter__.return_value = MagicMock()
    
    with patch("Bio.Entrez.read") as mock_read:
        # Pierwsze wywołanie to wynik esearch, drugie to efetch
        mock_read.side_effect = [sample_esearch_result, sample_efetch_result]
        
        df = fetch_pubmed_data(2023, mock_config)

    assert not df.empty
    assert len(df) == 2
    
    row1 = df.iloc[0]
    assert row1["pmid"] == "12345"
    assert row1["title"] == "CRISPR advances"
    assert row1["journal"] == "Nature"
    assert row1["doi"] == "10.1038/nature123"
    
    row2 = df.iloc[1]
    assert row2["pmid"] == "67890"
    assert row2["doi"] is None


@patch("Bio.Entrez.esearch")
def test_no_records_found(mock_esearch, mock_config):
    with patch("Bio.Entrez.read") as mock_read:
        mock_read.return_value = {"IdList": []}
        
        df = fetch_pubmed_data(2023, mock_config)
        
    assert isinstance(df, pd.DataFrame)
    assert df.empty


@patch("Bio.Entrez.efetch")
@patch("Bio.Entrez.esearch")
def test_determinism(mock_esearch, mock_efetch, mock_config, sample_esearch_result, sample_efetch_result):
    with patch("Bio.Entrez.read") as mock_read:
        # Symulujemy te same dane dla dwóch uruchomień
        mock_read.side_effect = [
            sample_esearch_result, sample_efetch_result,
            sample_esearch_result, sample_efetch_result
        ]
        
        df1 = fetch_pubmed_data(2023, mock_config)
        df2 = fetch_pubmed_data(2023, mock_config)
        
        pd.testing.assert_frame_equal(df1, df2)