# MCL
# Language: R
# Input: CSV (network)
# Output: prefix 
# Tested with: PluMA 1.1, R 4.0.0

PluMA plugin to take a weighted network and run Markov Clustering (van Dongen, 2000).
The plugin accepts as input a CSV file representing the network, with nodes represented
as both rows and columns and entry (i, j) is the weight of the edge from node i to node j.
The plugin then outputs clusters as separate CSV files using the output prefix, in the following format:

"","x"
"1","Family.Porphyromonadaceae.0005"
"2","Family.Porphyromonadaceae.0006"
"3","Family.Lachnospiraceae.0045"
"","x"
"1","Family.Lachnospiraceae.0001"
"2","Family.Lachnospiraceae.0002"
"3","Family.Lachnospiraceae.0007"
"4","Family.Lachnospiraceae.0013"
"5","Family.Ruminococcaceae.0003"
"6","Family.Lachnospiraceae.0018"
"7","Family.Ruminococcaceae.0006"
"8","Family.Lachnospiraceae.0044"
"9","Order.Clostridiales.0004"
"10","Family.Lachnospiraceae.0057"
"11","Dorea.0001"
"12","Oscillibacter.0004"

Note each cluster is separated by the delimiter "","x".

Updated August 2019: The plugin has been extended for first- and second- level clustering
First-level clustering files end with *.1.csv and second-level with *.2.csv
Individual cluster data is also output with the second-level clustering
