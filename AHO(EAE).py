#!/usr/bin/env python
# coding: utf-8

# In[6]:


from Bio import Entrez
from ahocorasick import Automaton
from collections import defaultdict

# Provide an email to Entrez in case they need to contact you
Entrez.email = "your_email_here@example.com"  # Change to your actual email

# Phrases to determine if it is pro_inflammatory or anti_inflammatory 
pro_inflammatory_phrases = [
    'pro-inflammatory', 'increases inflammation', 'induces inflammation',
    'elevates cytokine', 'upregulates NF-kB', 'activates macrophages',
    'stimulates immune response', 'enhances immune reaction'
]

anti_inflammatory_phrases = [
    'anti-inflammatory', 'reduces inflammation', 'suppresses inflammation',
    'decreases cytokine production', 'downregulates NF-kB', 'inhibits macrophage activation',
    'dampens immune response', 'attenuates immune reaction'
]

def fetch_article_ids(search_terms, num_results=10000):
    combined_search_term = " OR ".join(['"%s"' % term for term in search_terms])
    handle = Entrez.esearch(db="pubmed", term=combined_search_term, retmax=num_results)
    result = Entrez.read(handle)
    handle.close()
    return result["IdList"]

def fetch_abstracts(article_ids):
    handle = Entrez.efetch(db="pubmed", id=','.join(article_ids), rettype="medline", retmode="xml")
    records = Entrez.read(handle)
    abstracts = [record['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [""])[0] for record in records['PubmedArticle']]
    handle.close()
    return abstracts

def infer_inflammation_role(abstract, gene_protein_name):
    # Split abstract into words
    words = abstract.split()
    # Find the position(s) of the gene/protein name in the abstract
    positions = [i for i, word in enumerate(words) if word == gene_protein_name]
    
    # Define the window size for checking phrases before and after the gene/protein name
    window_size = 5
    
    # Initialize counts
    pro_count, anti_count = 0, 0
    
    for position in positions:
        # Check the surrounding words within the window size
        start = max(0, position - window_size)
        end = min(len(words), position + window_size + 1)
        surrounding = " ".join(words[start:end])
        
        # Count occurrences of pro and anti-inflammatory phrases
        pro_count += sum(phrase in surrounding for phrase in pro_inflammatory_phrases)
        anti_count += sum(phrase in surrounding for phrase in anti_inflammatory_phrases)
    
    # Infer the role based on the higher count
    if pro_count > anti_count:
        return 'pro-inflammatory'
    elif anti_count > pro_count:
        return 'anti-inflammatory'
    else:
        return 'unknown'

def build_aho_corasick_tree(keywords):
    A = Automaton()
    for keyword in keywords:
        A.add_word(keyword, keyword)
    A.make_automaton()
    return A

def search_abstracts(abstracts, automaton):
    found_keywords = defaultdict(list)
    for idx, abstract in enumerate(abstracts):
        for end_index, found in automaton.iter(abstract):
            keyword = found
            found_keywords[keyword].append(idx + 1)  # Use idx+1 to match abstract numbering
    return found_keywords

# Change the search query
search_terms = ["EAE", "experimental autoimmune encephalitis", "experimental autoimmune encephalomyelitis"]
print(f"Querying PubMed for articles related to {', '.join(search_terms)}...")
article_ids = fetch_article_ids(search_terms)
print(f"Fetching abstracts for {len(article_ids)} articles...")
abstracts = fetch_abstracts(article_ids)

print("\nLoading all genes and proteins from the file...")
genes_proteins = []
with open("/Users/sanjanakavula/Downloads/mart_export.txt", "r") as f:
    next(f)  # Skip the header line
    for line in f:
        _, gene_protein_name = line.strip().split(',')
        if gene_protein_name and len(gene_protein_name) > 2:  # Exclude empty or very short names
            genes_proteins.append(gene_protein_name)
print(f"Loaded {len(genes_proteins)} genes and proteins.")

print("\nBuilding the Aho-Corasick tree...")
automaton = build_aho_corasick_tree(genes_proteins)

print("\nSearching abstracts for gene and protein mentions using Aho-Corasick...")
found_keywords_in_abstracts = search_abstracts(abstracts, automaton)

# Frequency analysis
gene_protein_frequency = defaultdict(int)
for keyword, abstract_nums in found_keywords_in_abstracts.items():
    gene_protein_frequency[keyword] += len(abstract_nums)

# Print the frequency analysis results
print("\nGene and Protein Frequency Analysis:")
sorted_by_frequency = sorted(gene_protein_frequency.items(), key=lambda x: x[1], reverse=True)
for name, freq in sorted_by_frequency:
    print(f"{name}: {freq} mentions")

#Print Key words in genes along with predicted role
for keyword, abstract_nums in found_keywords_in_abstracts.items():
    if abstract_nums:  # Print only if there are abstract numbers
        inferred_roles = [infer_inflammation_role(abstracts[num-1], keyword) for num in abstract_nums]
        print(f"'{keyword}' is mentioned in abstracts: {abstract_nums}")
        # Count the roles for the current keyword
        role_counts = defaultdict(int)
        for role in inferred_roles:
            role_counts[role] += 1
        # Print the role counts for the current keyword
        for role, count in role_counts.items():
            print(f"  Inferred role - {role}: {count} times")


# In[ ]:




