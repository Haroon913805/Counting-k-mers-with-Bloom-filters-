from pybloom_live import BloomFilter
import os
import tarfile

# a) Create Bloom Filters for Viral Genomes
def create_bloom_filters(genomes_path):
    print("Creating Bloom filters for viral genomes...")
    viral_genome_bloom_filters = {}
    try:
        with tarfile.open(genomes_path, "r:gz") as tar:
            for member in tar.getmembers():
                if member.isreg():  # Check if member is a file
                    genome_name = os.path.basename(member.name).split('.')[0]
                    print(f"Processing genome: {genome_name}")
                    with tar.extractfile(member) as genome_file:
                        genome_sequence = genome_file.read()  # Read in binary mode
                        k_mers = [genome_sequence[i:i+10] for i in range(len(genome_sequence)-9) if b'N' not in genome_sequence[i:i+10]]
                        bloom_filter = BloomFilter(capacity=len(k_mers), error_rate=0.01)
                        for k_mer in k_mers:
                            bloom_filter.add(k_mer)  # Keep k_mers as binary
                        viral_genome_bloom_filters[genome_name] = bloom_filter
        print("Bloom filters creation completed.")
    except FileNotFoundError:
        print(f"Error: File '{genomes_path}' not found.")
    except tarfile.ReadError:
        print(f"Error: Failed to open or read file '{genomes_path}'.")
    except Exception as e:
        print(f"Error: An unexpected error occurred: {e}")

    return viral_genome_bloom_filters

# b) Assign Reads to Viral Genomes
def assign_reads_to_genomes(reads_path, viral_genome_bloom_filters):
    print("Assigning reads to viral genomes...")
    try:
        with open(reads_path, "r") as reads_file:
            for line in reads_file:
                if line.startswith(">"):  # Skip header lines
                    continue
                read_sequence = line.strip()
                read_k_mers = [read_sequence[i:i+10] for i in range(len(read_sequence)-9) if 'N' not in read_sequence[i:i+10]]
                assigned_genome = None
                max_matches = 0
                for genome_name, bloom_filter in viral_genome_bloom_filters.items():
                    matches = sum(bloom_filter.check(k_mer.encode('utf-8')) for k_mer in read_k_mers)  # Encode k_mer to binary
                    if matches > max_matches:
                        assigned_genome = genome_name
                        max_matches = matches
                if assigned_genome:
                    print(f"Read assigned to {assigned_genome}")
        print("Read assignment completed.")
    except FileNotFoundError:
        print(f"Error: File '{reads_path}' not found.")
    except Exception as e:
        print(f"Error: An unexpected error occurred: {e}")

# Main code
if __name__ == "__main__":
    genomes_path = "genomes .tar.gz"
    reads_path = "reads.fa"

    # a) Create Bloom Filters for Viral Genomes
    viral_genome_bloom_filters = create_bloom_filters(genomes_path)

    # b) Assign Reads to Viral Genomes
    assign_reads_to_genomes(reads_path, viral_genome_bloom_filters)