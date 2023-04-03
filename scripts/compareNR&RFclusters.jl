using PhyloNetworks, PhyloPlots
include("RFutilities.jl")
using RCall

## EXTRACT SOFTWIRED CLUSTERS OF INFERRED NETRAX NETWORKS ##
l2r_str = "((((BHV5:0.00899701)#6:0.0069709::0.939003)#8:0.0069709::0.898668)#5:0.315321::0.927192,(((((((((B589:0.000791084)#3:0.00157673::0.708934,(K22:0.000616704)#4:0.00117436::0.852769)38:0.000813901,SP1777:0.00203986)57:0.00128277,((((((Cooper:3.34627e-05)#7:7.9853e-05::0.163911)#0:0.00195626::0.901044,(216_II:1e-06)#2:0.000444769::0.801882)99:0.00120403,((Titanium_IBR_MLV_vaccine:8.3631e-05,((C46:5.74592e-06,(#0:1e-06::0.0989555,#8:1e-06::0.101332):0.500022)34:1.04782e-05,#7:2e-06::0.836089):7.85011e-05):0.000147567,C33:0.000202677):0.000904741)52:0.00200274)#1:0.000525991::0.797331,#6:1e-06::0.0609966):0.00157797)36:0.0016303,#3:0.00142325::0.291066):0.00338619,#4:0.000273377::0.147231):0.0107843,#2:0.00496123::0.198118):0.00229439,#1:0.0026344::0.202669):0.00645918,#5:1e-06::0.0728081):1e-06)Root;";
# Pre-process newick string so that hybrid node names have the form #H<name>
# E.g. #H1 instead of #1
l2r_net = replace(l2r_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
    readTopology; # l2r_net.numHybrids: 9
l2r_9ret = deepcopy(l2r_net);
shrink3cycles!(l2r_9ret);
# length(displayedTrees(l2r_9ret,-1.0)) # 144
# log(2,144) ≈ 7.2

r2l_str = "(((216_II:0.00105963)#0:1e-06::0.0106441,(BHV5:0.00388591)#1:0.17141::0.959828):0.198097,((((((Titanium_IBR_MLV_vaccine:8.61055e-05,(C46:2.10833e-05,(Cooper:0.000103757)#7:1e-06::0.889918):8.74804e-05)48:0.000160024,C33:0.000208295):1e-06)#3:0.0042661::0.17765,((((#0:0.000718469::0.989356,#3:0.000965679::0.82235)51:0.00208014,#1:1e-06::0.0401716):9.20093e-05)#2:1e-06::0.419125,#7:0.00478395::0.110082)85:0.00283014)42:0.00283014,(((SP1777:0.00174693,(((#2:0.000429798::0.580875,(K22:0.00162881)#4:0.000450355::0.322596):0.000653624,(B589:0.00113845)#6:0.00236132::0.469234)96:0.000479139)#5:0.000620419::0.827606)83:0.000860251,#4:1e-06::0.677404):0.000309128,#6:1e-06::0.530766):0.00340796):0.00667725,#5:1e-06::0.172394):0.00211855)Root;";
r2l_net = replace(r2l_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
    readTopology; # r2l_net.numHybrids: 8
r2l_8ret = deepcopy(r2l_net);
shrink3cycles!(r2l_8ret);
# length(displayedTrees(r2l_8ret,-1.0)) # 224
# log(2,224) ≈ 7.8

taxon_labels = tipLabels(l2r_9ret);

clusters_l2r = foldl(union, # union of softwired clusters (encoded as integers)
    map(t -> clust2int(hardwiredClusters(t,taxon_labels))[1],
    displayedTrees(l2r_9ret,0.0));init=Int[]) |> sort
# [6, 7, 18, 22, 23, 24, 31, 55, 63, 96, 119, 127, 160, 183, 191, 192, 224, 232, 247, 255, 258, 262, 263, 272, 274, 278, 279, 287, 288, 311, 319, 352, 375, 383, 416, 439, 447, 480, 488, 503]
# length(clusters_l2r): 40


clusters_r2l = foldl(union,
    map(t -> clust2int(hardwiredClusters(t,taxon_labels))[1],
    displayedTrees(r2l_8ret,0.0));init=Int[]) |> sort
# [6, 7, 15, 18, 22, 23, 24, 31, 39, 40, 47, 48, 55, 63, 71, 72, 79, 87, 95, 96, 103, 104, 111, 112, 119, 127, 135, 136, 143, 151, 159, 160, 167, 168, 175, 176, 183, 191, 192, 199, 200, 207, 215, 223, 224, 231, 232, 239, 240, 247, 255, 263, 264, 271, 272, 279, 280, 287, 288, 295, 296, 303, 311, 319, 320, 327, 328, 335, 343, 351, 352, 359, 360, 367, 375, 383, 384, 391, 392, 399, 407, 415, 416, 423, 424, 431, 439, 447, 448, 455, 456, 463, 471, 479, 480, 487, 488, 495, 503]
# length(clusters_r2l): 99

################################################################################

## EXTRACT SOFTWIRED CLUSTERS OF INFERRED RF-NET NETWORKS ##
l2r_7ret_rf = readMultiTopology("../results/bhv/rfnet/97genes_blksize1500_set1c_ogrooted.newick")[8];
# displayedTrees(l2r_7ret_rf,-1.0) |> length # 128
# log(2,128) # 7

r2l_7ret_rf = readMultiTopology("../results/bhv/rfnet/97genes_blksize1500_rev_set1c_ogrooted.newick")[8];
# displayedTrees(r2l_7ret_rfnet,-1.0) |> length # 120
# log(2,120) ≈ 6.9

clusters_l2r_rf = foldl(union,
    map(t -> clust2int(hardwiredClusters(t,taxon_labels))[1],
    displayedTrees(l2r_7ret_rf,-1.0));init=Int[]) |> sort
# [3, 5, 6, 7, 17, 18, 19, 20, 21, 22, 23, 31, 63, 95, 96, 127, 160, 191, 192, 223, 224, 255]
# length(clusters_l2r_rf) # 22

clusters_r2l_rf = foldl(union,
    map(t -> clust2int(hardwiredClusters(t,taxon_labels))[1],
    displayedTrees(r2l_7ret_rf,-1.0));init=Int[]) |> sort
# [3, 5, 6, 7, 17, 18, 19, 20, 21, 22, 23, 31, 63, 96, 127, 160, 191, 192, 224, 255]
# length(clusters_r2l_rf) # 20

################################################################################

## COMMON SOFTWIRED CLUSTERS
intersect(clusters_l2r,clusters_l2r_rf)
# 6, 7, 18, 22, 23, 31, 63, 96, 127, 160, 191, 192, 224, 255]

intersect(clusters_r2l,clusters_r2l_rf)
# [6, 7, 18, 22, 23, 31, 63, 96, 127, 160, 191, 192, 224, 255]

for c in intersect(clusters_l2r,clusters_l2r_rf)
    printCluster(c,9,taxon_labels)
end
# ["Titanium_IBR_MLV_vaccine", "C46"]
# ["Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["Cooper", "C46"]
# ["Cooper", "Titanium_IBR_MLV_vaccine", "C46"]
# ["Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["SP1777", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["K22", "SP1777"]
# ["K22", "SP1777", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["B589", "SP1777"]
# ["B589", "SP1777", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["B589", "K22"]
# ["B589", "K22", "SP1777"]
# ["B589", "K22", "SP1777", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]

## SOFTWIRED CLUSTERS NOT SHARED
for c in symdiff(clusters_l2r,clusters_l2r_rf)
    printCluster(c,9,taxon_labels)
end
# symdiff(clusters_l2r,clusters_l2r_rf) |> length # 34
# union(clusters_l2r,clusters_l2r_rf) |> length # 48
# ["Cooper", "216_II"]
# ["SP1777", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["K22", "SP1777", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["B589", "SP1777", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["B589", "K22", "SP1777", "216_II"]
# ["B589", "K22", "SP1777", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "C46"]
# ["BHV5", "Titanium_IBR_MLV_vaccine", "C46"]
# ["BHV5", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "Cooper"]
# ["BHV5", "Cooper", "C46"]
# ["BHV5", "Cooper", "Titanium_IBR_MLV_vaccine", "C46"]
# ["BHV5", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "SP1777"]
# ["BHV5", "SP1777", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "SP1777", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "K22", "SP1777"]
# ["BHV5", "K22", "SP1777", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "K22", "SP1777", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "B589", "SP1777"]
# ["BHV5", "B589", "SP1777", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "B589", "SP1777", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["BHV5", "B589", "K22", "SP1777"]
# ["BHV5", "B589", "K22", "SP1777", "216_II"]
# ["BHV5", "B589", "K22", "SP1777", "Cooper", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["C46", "C33"]
# ["Titanium_IBR_MLV_vaccine", "C33"]
# ["Cooper", "C33"]
# ["Cooper", "C46", "C33"]
# ["Cooper", "Titanium_IBR_MLV_vaccine"]
# ["Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["K22", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]
# ["B589", "K22", "Cooper", "216_II", "Titanium_IBR_MLV_vaccine", "C46", "C33"]

for c in symdiff(clusters_r2l,clusters_r2l_rf)
    printCluster(c,9,taxon_labels)
end
# symdiff(clusters_r2l,clusters_r2l_rf) |> length # 91
# union(clusters_r2l,clusters_r2l_rf) |> length # 105