
from pybloom_live import BloomFilter
import random


n = 1000000  
error_rate = 0.01  


bloom_filter = BloomFilter(capacity=n, error_rate=error_rate)

def generate_random_kmer(k):
    return ''.join(random.choices('ACGT', k=k))


print("Adding k-mers to the Bloom filter...")
for i in range(n):
    if i % 100000 == 0:
        print(f"Added {i} k-mers")
    kmer = generate_random_kmer(10)  
    bloom_filter.add(kmer)

print("Finished adding k-mers.")

def count_frequent_kmers(bloom_filter, n):
    kmer_counts = {}
    frequent_kmers_count = 0
    
    print("Counting k-mers...")
    for i in range(n):
        if i % 100000 == 0:
            print(f"Processed {i} k-mers")
        kmer = generate_random_kmer(10)
