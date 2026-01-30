import argparse
import pandas as pd
from plotnine import *
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser(
    description=(
        'Process FASTQ files to compute read-length and 5′ nucleotide statistics. '
        'Outputs include per-replicate summaries and an averaged summary in TSV format.'
    )
)


parser.add_argument(
    '-l', '--length-range',
    type=int,
    nargs=2,
    metavar=('MIN', 'MAX'),
    help=(
        'Minimum and maximum read length (inclusive) to include in the analysis. '
        'Example: -l 18 40'
    )
)

parser.add_argument(
    '-p', '--prefix',
    required=True,
    help=(
        'Sample prefix used to construct FASTQ filenames. Files must be named like PREFIX_1.fastq, PREFIX_2.fastq, ...'
    )
)

parser.add_argument(
    '-r', '--rep-number',
    type=int,
    default=4,
    help='Number of replicates per sample (default: 4)'
)

args = parser.parse_args()

lower, upper = args.length_range
replicate_prefix = args.prefix
replicate_number = args.rep_number

# Process the replicate prefix
n = (replicate_number + 1)
replicate_files = [f"{replicate_prefix}_{i}.fastq" for i in range(1, n)]
#####


#Create a list of fastq file names
fastq_files = replicate_files

# Create an empty dictionary of dictionaries to store read length statistics for each replicate
read_length_stats = defaultdict(lambda: defaultdict(int))

# Loop over fastq files
for fastq_file in fastq_files:
    # Process the fastq file
    with open(fastq_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fastq'):

            # Get length and starting nucleotide of read
            starting_nucleotide = record.seq[0]
            read_length = len(record.seq)

            # Increment count of reads with this length and starting nucleotide for this replicate
            read_length_stats[fastq_file][(starting_nucleotide, read_length)] += 1

for replicate, counts in read_length_stats.items():
    df = pd.DataFrame([(key[0], key[1], value) for key, value in counts.items()], columns=['Starting Nucleotide', 'Length', 'Raw Count'])
    filt_size = (df['Length'] >= lower) & (df['Length'] <= upper)
    df = df.loc[filt_size]
     # Normalize reads based on sequencing depth
    df['Reads (TPM)'] = ((df['Raw Count'] / df['Raw Count'].sum()) * 1000000)
    tsv_file = replicate.replace('.fastq', '.tsv')
    df.to_csv(tsv_file, sep='\t', index=False)

# Combine the read length statistics for all replicates
result = defaultdict(int)
for replicate_data in read_length_stats.values():
    for key, value in replicate_data.items():
        result[key] += value

# Calculate average read length statistics across all replicates
result = defaultdict(int)
for replicate_data in read_length_stats.values():
    for key, value in replicate_data.items():
        result[key] += value
result = {key: round(value / len(read_length_stats)) for key, value in result.items()}

# Convert the read length statistics to a pandas DataFrame
data = [(key[0], key[1], value) for key, value in result.items()]
df = pd.DataFrame(data, columns=['Starting_Nucleotide', 'Length', 'Average Reads'])

# Normalize reads based on sequencing depth
df['Average Reads (TPM)'] = ((df['Average Reads'] / df['Average Reads'].sum()) * 1000000)

# Filter size distribution based on your interests
filt_size = (df['Length'] >= lower) & (df['Length'] <= upper)
df = df.loc[filt_size]

# Write to a csv
output_filename = replicate_prefix
df.to_csv(output_filename + 'AVG.tsv', sep="\t")

# plot data
plot = (ggplot(df, aes('Length', 'Average Reads (TPM)', fill='Starting_Nucleotide'))
+ geom_col()
+ theme_classic()
#scale_fill_manual(values=['#D5CEC3','#F0C988', '#B5C48B', '#E89C82'])
+ labs(title="Average Size and 5' Nucleotide Distribution", x="Read Length (nt)", y='Reads per million (RPM)')
+ theme(plot_title=element_text(size=16, face='bold'),  # Change title appearance
            axis_text=element_text(size=12),  # Change axis text appearance
            axis_title=element_text(size=14, face='bold'))  # Change axis title appearance
)


plot
plot.save(output_filename + str(lower) + str(upper) + 'AVG.pdf', height=6, width=8)
