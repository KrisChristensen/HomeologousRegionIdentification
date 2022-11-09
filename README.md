# HomeologousRegionIdentification
A pipeline to identify homologous or homeologous regions within or between genomes.  A masked genome is aligned to itself or another genome using blast and then the alignments are filtered.  The output is in a format suitable for generating circos plots.

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- requirements -->
## Requirements

These scripts have been tested with Python 3.
The scripts require the following programs and files.

Programs:<br /><br />
&nbsp;&nbsp;&nbsp;blast (must be in format #6, see example below)<br />
    
Files:<br /><br />
&nbsp;&nbsp;&nbsp;A masked genome file or two depending on use<br />

<!-- usage -->
## Usage

1) Generate hard-masked genome (from NCBI, lowercase basepairs are softmasked):<br /><br />

&nbsp;&nbsp;&nbsp;&nbsp;example:
&nbsp;&nbsp;&nbsp;sed -e '/^>/! s/[[:lower:]]/N/g' GCF_002021735.2_Okis_V2_genomic.fna > GCF_002021735.2_Okis_V2_genomic.masked.fna<br /><br />
    
2) Generate blast database:<br /><br />

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;example:
&nbsp;&nbsp;&nbsp;&nbsp;makeblastdb -in GCF_002021735.2_Okis_V2_genomic.masked.fna -dbtype nucl<br /><br />
    
3) Align blast database to self:<br /><br />

&nbsp;&nbsp;&nbsp;&nbsp;example:
&nbsp;&nbsp;&nbsp;&nbsp;blastn -task megablast -db GCF_002021735.2_Okis_V2_genomic.masked.fna -query GCF_002021735.2_Okis_V2_genomic.masked.fna -out GCF_002021735.2.vs.self.aln -outfmt 6 -perc_identity 80 -max_hsps 40000<br /><br />
    
4) Identify homeologous regions in genome:<br /><br />
&nbsp;&nbsp;&nbsp;python General_linear_filter_fmt6.v1.3.py -aln GCF_002021735.2.vs.self.aln -gap 100000 -min 10000 -print no 2> Gap100K.Min10K.txt<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python General_linear_filter_fmt6.v1.3.py -h<br />
&nbsp;&nbsp;&nbsp;note: requires Linear_Alignments_v4.py and GeneralOverlap_v1.py in same working directory
    
5) Output in Circos plot format and filter small homeologous regions:<br /><br />
&nbsp;&nbsp;&nbsp;python CircosOutput.v1.1.py -input Gap100K.Min10K.txt -tbl ChrNameColor.txt -tLen 100000 > Gap100K.Min10K.circos.txt 2> Gap100K.Min10K.pid.circos.txt<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python CircosOutput.v1.1.py -h<br />
&nbsp;&nbsp;&nbsp;The tbl is a tab-delimited file with the chromosome name in the first column. New chromosome in the second column (can be the same as first column) and the color for the circos plot. See example file in repository.
    

<!-- license -->
## License 

Distributed under the MIT License.
