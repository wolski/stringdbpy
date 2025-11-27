# STRING-DB ORA (Over-Representation Analysis) API Documentation

This document describes how we call the STRING-DB API for functional enrichment analysis with a custom background.

## Overview

We use three STRING-DB API endpoints:
1. **`/json/get_string_ids`** - Map protein/gene identifiers to STRING IDs
2. **`/json/enrichment`** - Run functional enrichment analysis
3. **`/json/get_link`** - Get a link to the STRING-DB network visualization

**Base URL**: `https://version-12-0.string-db.org/api`

---

## Step 1: Map Identifiers to STRING IDs

Before calling the enrichment endpoint with a custom background, identifiers must be mapped to STRING IDs (format: `{taxon}.{ENSEMBL_protein_id}`).

### Endpoint
```
GET /json/get_string_ids
```

### Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `identifiers` | string | Protein/gene identifiers separated by `%0D` (carriage return) |
| `species` | int | NCBI Taxon ID (e.g., 10090 for mouse, 9606 for human) |
| `caller_identity` | string | Your application identifier |

### Example Request
```bash
curl "https://version-12-0.string-db.org/api/json/get_string_ids?identifiers=Brca1%0DTp53&species=10090"
```

### Example Response
```json
[
  {
    "queryIndex": 0,
    "queryItem": "Brca1",
    "stringId": "10090.ENSMUSP00000017290",
    "ncbiTaxonId": 10090,
    "taxonName": "Mus musculus",
    "preferredName": "Brca1",
    "annotation": "Breast cancer type 1 susceptibility protein..."
  },
  {
    "queryIndex": 1,
    "queryItem": "Tp53",
    "stringId": "10090.ENSMUSP00000104298",
    "ncbiTaxonId": 10090,
    "taxonName": "Mus musculus",
    "preferredName": "Trp53",
    "annotation": "Cellular tumor antigen p53..."
  }
]
```

---

## Step 2: Run Enrichment Analysis

### Endpoint
```
POST /json/enrichment
```

**Note**: We use POST to avoid URL length limits when sending many identifiers.

### Parameters
| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `identifiers` | string | Yes | STRING IDs of significant proteins, separated by `%0D` |
| `background_string_identifiers` | string | No | STRING IDs for custom background, separated by `%0D` |
| `species` | int | Yes | NCBI Taxon ID |
| `caller_identity` | string | Yes | Your application identifier |

### Example Request (with custom background)
```bash
curl -X POST "https://version-12-0.string-db.org/api/json/enrichment" \
  -d "identifiers=10090.ENSMUSP00000017290%0D10090.ENSMUSP00000104298" \
  -d "background_string_identifiers=10090.ENSMUSP00000017290%0D10090.ENSMUSP00000104298%0D10090.ENSMUSP00000000001" \
  -d "species=10090" \
  -d "caller_identity=my_application"
```

### Example Request (without custom background - uses genome-wide)
```bash
curl "https://version-12-0.string-db.org/api/json/enrichment?identifiers=Brca1%0DTp53&species=10090"
```

### Example Response
```json
[
  {
    "category": "COMPARTMENTS",
    "term": "GOCC:0005654",
    "number_of_genes": 92,
    "number_of_genes_in_background": 586,
    "ncbiTaxonId": 10090,
    "inputGenes": ["10090.ENSMUSP00000001258", "10090.ENSMUSP00000001326", ...],
    "preferredNames": ["Uhrf1", "Sp1", ...],
    "p_value": 2.36e-8,
    "fdr": 0.0000304,
    "description": "Nucleoplasm"
  },
  ...
]
```

### Response Fields
| Field | Description |
|-------|-------------|
| `category` | Enrichment category (COMPARTMENTS, Component, Process, Function, KEGG, RCTM, WikiPathways, DISEASES, TISSUES) |
| `term` | Term identifier (e.g., GO:0005654, KEGG:hsa04110) |
| `number_of_genes` | Number of input genes in this term |
| `number_of_genes_in_background` | Number of background genes in this term |
| `ncbiTaxonId` | Species taxon ID |
| `inputGenes` | List of STRING IDs that matched this term |
| `preferredNames` | List of gene names |
| `p_value` | Enrichment p-value |
| `fdr` | False discovery rate (Benjamini-Hochberg corrected) |
| `description` | Human-readable term description |

---

## Step 3: Get Network Link (Optional)

### Endpoint
```
POST /json/get_link
```

**Note**: Use POST for many identifiers to avoid URL length limits.

### Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `identifiers` | string | STRING IDs separated by `%0D` |
| `species` | int | NCBI Taxon ID |
| `caller_identity` | string | Your application identifier |

### Example Request
```bash
curl -X POST "https://version-12-0.string-db.org/api/json/get_link" \
  -d "identifiers=10090.ENSMUSP00000017290%0D10090.ENSMUSP00000104298" \
  -d "species=10090" \
  -d "caller_identity=my_application"
```

### Example Response
```json
["https://version-12-0.string-db.org/cgi/link?to=94087B98627657B6"]
```

---

## Feature Request: Permanent Link to ORA Results

Currently, the `/json/get_link` endpoint returns a link to the **network visualization**, not to the **enrichment/ORA results page**.

**Requested feature**: An API endpoint or parameter that returns a permanent/shareable link to the ORA enrichment results page (similar to what users see when they run enrichment analysis through the web interface).

### Current Workflow
1. Map identifiers → STRING IDs
2. Call `/json/enrichment` with `background_string_identifiers`
3. Receive JSON results directly (no URL to results page)
4. `/json/get_link` only provides network visualization link

### Desired Workflow
1. Map identifiers → STRING IDs
2. Call enrichment endpoint
3. Receive both:
   - JSON results
   - Permanent URL to view enrichment results on STRING-DB website

---

## Python Implementation

See `src/string_gsea/scripts/string_ora_run.py` for the full implementation.

### Key Code Snippets

#### Mapping Identifiers
```python
url = "https://version-12-0.string-db.org/api/json/get_string_ids"
params = {
    "identifiers": "\r".join(identifiers),  # \r = %0D when URL encoded
    "species": species,
    "caller_identity": caller_identity,
}
resp = requests.get(url, params=params)
data = resp.json()
# Extract: {query_item: string_id}
```

#### Running Enrichment with Background
```python
url = "https://version-12-0.string-db.org/api/json/enrichment"
params = {
    "identifiers": "\r".join(significant_string_ids),
    "background_string_identifiers": "\r".join(background_string_ids),
    "species": species,
    "caller_identity": caller_identity,
}
resp = requests.post(url, data=params)  # POST for large payloads
results = resp.json()
```

---

## References

- STRING-DB API Documentation: https://string-db.org/cgi/help.pl?subpage=api
- STRING-DB Enrichment: https://string-db.org/cgi/help.pl?subpage=api#getting-functional-enrichment
