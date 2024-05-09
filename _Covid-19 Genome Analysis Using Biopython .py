#!/usr/bin/env python
# coding: utf-8

# In[15]:


from Bio import Entrez, SeqIO
Entrez.email = "" 
handle = Entrez.efetch(db="nucleotide", id="MK967708.1", rettype="gb", retmode="text")
recs = list(SeqIO.parse(handle, 'gb'))
handle.close()


# In[16]:


recs


# In[17]:


covid_rna = recs[0].seq


# In[18]:


covid_rna


# In[19]:


print(f'The genome of Covid-19 consists of {len(covid_rna)} nucleotides.')


# In[20]:


#moleclar weight
from Bio.SeqUtils import molecular_weight
molecular_weight(covid_rna)


# In[21]:


# GC content - higher GC content implies more stable molecule due to G and C forming triple hydrogen bonds
from Bio.SeqUtils import gc_fraction
gc_fraction(covid_rna)
 


# In[22]:


#Distribution of neuclotides


# In[23]:


count_nucleotides = {
    'A': covid_rna.count('A'),
    'T': covid_rna.count('T'),
    'C': covid_rna.count('C'),
    'G': covid_rna.count('G')
}


# In[24]:


print(count_nucleotides)


# In[25]:


import matplotlib.pyplot as plt
width = 0.5
plt.bar(count_nucleotides.keys(), count_nucleotides.values(), width, color=['b', 'r', 'g', 'c'])
plt.xlabel('Nucleotide')
plt.ylabel('Frequency')
plt.title('Nucleotide Frequency')


# In[26]:


#TRANSCRIPTION


# In[27]:


covid_mrna = covid_rna.transcribe()
covid_mrna


# In[28]:


#TRANSLATION


# In[29]:


covid_aa = covid_mrna.translate()
covid_aa


# In[30]:


#MOST COMMON AMINO ACIDS


# In[31]:


from collections import Counter
common_amino = Counter(covid_aa)
common_amino.most_common(10)


# In[32]:


del common_amino['*']

width = 0.5
plt.bar(common_amino.keys(), common_amino.values(), width, color=['b', 'r', 'm', 'c'])
plt.xlabel('Amino Acid')
plt.ylabel('Frequency')
plt.title('Protein Sequence Frequency')


# In[33]:


print(f"Covid-19's genome has {sum(common_amino.values())} amino acids")


# In[34]:


proteins = covid_aa.split('*')


# In[35]:


proteins[:5]


# In[37]:


print(f'We have {len(proteins)} proteins in the covid-19 genome')


# In[ ]:


#It's worth to mention that not all the amino acids sequences are proteins. Only the sequences with more than 20 amino acids code for functional proteins. The short amino acid sequences are oligopeptides and have other functionalities. Here, we will focus on the chains with more than 20 amino acid chains: Proteins.


# In[38]:


for protein in proteins:
    if len(protein) < 20:
        proteins.remove(protein)


# In[39]:


print(f'We have {len(proteins)} proteins with more than 20 amino acids in the covid-19 genome')


# In[40]:


top_5_proteins = sorted(proteins, key = len)


# In[41]:


top_5_proteins[-1]


# In[42]:


len(top_5_proteins[-1])


# In[43]:


with open("protein_seq.fasta", "w") as file:
    file.write(f">covid protein\n{top_5_proteins[-1]}")


# In[ ]:




