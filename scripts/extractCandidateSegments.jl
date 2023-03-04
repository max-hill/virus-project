using DelimitedFiles
using PhyloNetworks

# Preserves metadata (e.g. gene embeddings) that PhyloNetworks ignores
net_string = readlines("../results/bhv/rfnet/97genes_blksize1500_set1c_ogrooted.newick")[8];
numgenes = 97;

# Regex expression to capture hybrid name and supporting genes
rx_rfnet = r"#(?<hybrid>H[\d]+)\[&genes={(?<genes>[\d,]+)}\]";

# Search for all matches of `rx_rfnet` and return an iterator over the matches
hname2genes = eachmatch(rx_rfnet,net_string) |> collect;

d_hedge2support = Dict(); # Maps minor hybrid edges to gene support

# Maps minor hybrid edges to supporting "segments": "segments" are sets of genes
# that form a contiguous block. E.g. A minor hybrid edge supported by genes 13,
# 14, 15, is supported by the segment (13,15).
d_hedge2seg = Dict();

# Find all ranges of consecutive numbers from a sorted (in ascending order) array
function extractRanges(arr::Vector{Int})
    prev = arr[1]
    n = length(arr)
    if n == 1
        return [(prev,prev)]
    end
    ranges = Tuple{Int,Int}[] # represent each range as a tuple
    m = 1 # length of range ≥ 1
    for i = 2:n
        curr = arr[i]
        if (curr-prev) == 1
            m += 1 # increment the "length" variable
        else
            # build range using the previous element and "length" variable
            push!(ranges,(prev-m+1,prev))
            m = 1
        end
        prev = curr
    end
    return ranges
end

for i in eachindex(hname2genes)
    hname = hname2genes[i].captures[1] # hybrid name
    genes = parse.(Int64,
        split(hname2genes[i].captures[2],",")) # (sorted) gene numbers
    support = length(genes)
    if support < ceil(numgenes/2) # if minor edge
        d_hedge2support[hname] = support
        d_hedge2seg[hname] = extractRanges(genes) # segments
    end
end

d_hedge2support # out of 97 genes
# "H22" => 42
# "H18" => 43
# "H33" => 33
# "H20" => 16
# "H31" => 25
# "H27" => 29
# "H25" => 39

# Count number of segments for comparison
map(k -> k => length(d_hedge2seg[k]),collect(keys(d_hedge2seg)))
# "H22" => 21
# "H18" => 25
# "H33" => 20
# "H20" => 11
# "H31" => 19
# "H27" => 22
# "H25" => 26

# Look at segments of length > 1
map(k -> k => filter(t -> t[2] > t[1],d_hedge2seg[k]),collect(keys(d_hedge2seg)))
# "H22" => [(1, 2), (10, 12), (16, 20), (65, 67), (69, 70), (74, 82), (85, 87)]
# "H18" => [(13, 15), (34, 35), (51, 57), (72, 76), (85, 86)]
# "H33" => [(12, 16), (31, 32), (53, 54), (71, 74), (78, 80), (94, 95)]
# "H20" => [(11, 12), (26, 27), (38, 39), (59, 60)]
# "H31" => [(19, 20), (31, 32), (63, 64), (72, 73)]
# "H27" => [(4, 5), (15, 16), (19, 20), (45, 46), (74, 75), (87, 88)]
# "H25" => [(9, 10), (19, 20), (24, 25), (33, 38), (62, 63), (90, 93)]

# Look at segment lengths > 1
map(k -> k => map(t -> t[2]-t[1]+1,filter(t -> t[2] > t[1],d_hedge2seg[k])),
    collect(keys(d_hedge2seg)))
# "H22" => [2, 3, 5, 3, 2, 9, 3]
# "H18" => [3, 2, 7, 5, 2]
# "H33" => [5, 2, 2, 4, 3, 2]
# "H20" => [2, 2, 2, 2]
# "H31" => [2, 2, 2, 2]
# "H27" => [2, 2, 2, 2, 2, 2]
# "H25" => [2, 2, 2, 6, 2, 4]

# Regex expression to capture gene number and site range wrt the MSA
rx_part = r" gene_(?<gene>[\d]+)=(?<start>[\d]+)-(?<end>[\d]+)";

# Maps gene number to site range wrt the MSA (by reading partition file)
d_gene2siterange = Dict(
    map(
        function(str)
            g, s, e = parse.(Int64,match(rx_part,str))
            return (g,(s,e))
        end,
        readdlm("../data/bhv/part/97genes_blksize1500.txt",',',String)[:,2])
)

# Map segments of length > 1 (i.e. recombinant segment candidates) to site ranges
# wrt the MSA
map(
    function(k)
        map(t -> (d_gene2siterange[t[1]][1],d_gene2siterange[t[2]][2]),
            filter(t -> t[2] > t[1],d_hedge2seg[k]))
    end,collect(keys(d_hedge2seg)))
# [(1, 3000), (13501, 18000), (22501, 30000), (96001, 100500), (102001, 105000), (109501, 123000), (126001, 130500)]
# [(18001, 22500), (49501, 52500), (75001, 85500), (106501, 114000), (126001, 129000)]
# [(16501, 24000), (45001, 48000), (78001, 81000), (105001, 111000), (115501, 120000), (139501, 142500)]
# [(15001, 18000), (37501, 40500), (55501, 58500), (87001, 90000)]
# [(27001, 30000), (45001, 48000), (93001, 96000), (106501, 109500)]
# [(4501, 7500), (21001, 24000), (27001, 30000), (66001, 69000), (109501, 112500), (129001, 132000)]
# [(12001, 15000), (27001, 30000), (34501, 37500), (48001, 57000), (91501, 94500), (133501, 139500)]

# From BootScan, the recombinant segments from K22 detected in 216-II are located
# around (110000, 120000) and (136000, 144000).
# Above, we detect the candidate segments: (106501, 114000), (109501, 112500),
# (109501, 123000), (115501, 120000), (139501, 142500)

# From BootScan, the recombinant segment from BHV5 detected in 216-II corresponds
# to (81071, 82601).
# Above, we detect the candidate segment (75001, 85500).