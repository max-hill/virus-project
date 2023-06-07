###
using Pkg # to use functions that manage packages
Pkg.add("PhyloNetworks")
Pkg.add("PhyloPlots")
Pkg.add("Plots")
using PhyloNetworks
using PhyloPlots
using Plots

# For plotting networks and trees side-by-side
# Pkg.build("RCall")
Pkg.add("RCall")
using RCall
# displayedTrees(net0,0.41) |> length
###
#_______________________________________________________________________________
#
# Define starting trees and inferred networks.
#_______________________________________________________________________________
#_______________________________________________________________________________
#
# Define networks
#_______________________________________________________________________________

# The following are the 15 networks obtained from the multistart-titaniuam run,
# and the 15 starting trees obtained by IQtree. Each network was obtained by
# running netrax with the corresponding start tree.
# 
# Note that network 12 (obtained on the 13th run) has the highest likelihood. We
# had to rename the hybrid nodes so that they start with a character not 0. To
# do this, I replaced '#' with '#H'

### START CODE BLOCK
starting_tree_0=readTopology("(BHV5:0.14597,(SP1777:0.00219,B589:0.00259)NODE_000000599.00:0.00254,(216_II:0.00207,(C33:0.00011,(Titanium:0.00000,(C46:0.00000,Cooper:0.00000)NODE_000000287.00:0.00011)NODE_000000172.00:0.00011)NODE_0000000:0.00095)NODE_0000003100.00:0.00321)NODE_0000004100.00:0.00100;"); starting_tree_1=readTopology("(BHV5:0.07977,(SP1777:0.00394,B589:0.00254)NODE_000000596.00:0.00195,(216_II:0.00107,(C33:0.00031,(Titanium:0.00000,(C46:0.00000,Cooper:0.00000)NODE_000000274.00:0.00010)NODE_000000194.00:0.00031)NODE_0000000:0.00059)NODE_000000399.00:0.00382)NODE_0000004100.00:0.00100;"); starting_tree_2=readTopology("(BHV5:0.05776,((SP1777:0.00195,B589:0.00216)NODE_0000005100.00:0.00289,(216_II:0.00062,(C33:0.00010,(Titanium:0.00000,(C46:0.00000,Cooper:0.00010)NODE_000000283.00:0.00010)NODE_000000169.00:0.00010)NODE_0000000:0.00061)NODE_0000003100.00:0.00230)NODE_0000004100.00:0.05408)NODE_0000006:0.00100;"); starting_tree_3=readTopology("(BHV5:0.09174,(SP1777:0.00234,B589:0.00123)NODE_000000587.00:0.00152,(216_II:0.00143,(C33:0.00020,(Cooper:0.00000,(C46:0.00000,Titanium:0.00000)NODE_000000234.00:0.00000)NODE_0000001100.00:0.00051)NODE_0000000:0.00083)NODE_0000003100.00:0.00320)NODE_0000004100.00:0.00100;");   starting_tree_4=readTopology("(BHV5:0.16109,(SP1777:0.00196,B589:0.00378)NODE_000000547.00:0.00219,(216_II:0.00180,(C33:0.00000,(Cooper:0.00000,(C46:0.00011,Titanium:0.00011)NODE_000000238.00:0.00000)NODE_000000145.00:0.00000)NODE_0000000:0.00095)NODE_0000003100.00:0.00441)NODE_0000004100.00:0.00100;");   starting_tree_5=readTopology("(BHV5:0.08174,(SP1777:0.00111,B589:0.00061)NODE_000000598.00:0.00103,(216_II:0.00112,(C33:0.00020,(Titanium:0.00010,(C46:0.00000,Cooper:0.00000)NODE_000000246.00:0.00000)NODE_000000162.00:0.00010)NODE_0000000:0.00041)NODE_000000387.00:0.00296)NODE_0000004100.00:0.00100;");   starting_tree_6=readTopology("(BHV5:0.12529,(SP1777:0.00127,B589:0.00171)NODE_000000580.00:0.00092,(216_II:0.00273,(Titanium:0.00010,(C33:0.00021,(C46:0.00000,Cooper:0.00000)NODE_000000189.00:0.00010)NODE_0000000:0.00000)NODE_000000261.00:0.00140)NODE_0000003100.00:0.00366)NODE_000000480.00:0.00100;");   starting_tree_7=readTopology("(BHV5:0.07659,(SP1777:0.00072,B589:0.00175)NODE_000000591.00:0.00078,(216_II:0.00185,(C33:0.00010,(Titanium:0.00000,(C46:0.00000,Cooper:0.00000)NODE_000000238.00:0.00000)NODE_000000196.00:0.00031)NODE_0000000:0.00061)NODE_000000399.00:0.00253)NODE_0000004100.00:0.00100;");   starting_tree_8=readTopology("(BHV5:0.03335,(216_II:0.00158,((SP1777:0.00081,B589:0.00173)NODE_0000004100.00:0.00395,(Cooper:0.00000,(C33:0.00020,(C46:0.00000,Titanium:0.00000)NODE_000000153.00:0.00000)NODE_0000000:0.00000)NODE_000000243.00:0.00177)NODE_0000003100.00:0.00897)NODE_0000005100.00:0.02240)NODE_0000006:0.00100;");   starting_tree_9=readTopology("(BHV5:0.07048,(SP1777:0.00172,B589:0.00122)NODE_000000596.00:0.00130,(216_II:0.00213,(C33:0.00020,(Titanium:0.00010,(C46:0.00000,Cooper:0.00015)NODE_000000247.00:0.00000)NODE_000000172.00:0.00015)NODE_0000000:0.00172)NODE_0000003100.00:0.00216)NODE_0000004100.00:0.00100;");   starting_tree_10=readTopology("(B589:0.00210,BHV5:0.17992,(SP1777:0.00233,(216_II:0.00291,(C33:0.00049,(Titanium:0.00000,(C46:0.00000,Cooper:0.00264)NODE_000000267.00:0.00011)NODE_000000158.00:0.00017)NODE_0000000:0.00127)NODE_0000003100.00:0.00678)NODE_0000004100.00:0.00219)NODE_000000591.00:0.00100;");   starting_tree_11=readTopology("(B589:0.00285,BHV5:0.30906,(SP1777:0.00520,(216_II:0.00808,(Cooper:0.00407,(C33:0.00013,(C46:0.00025,Titanium:0.00049)NODE_000000158.00:0.00012)NODE_0000000:0.00026)NODE_000000268.00:0.00652)NODE_0000003100.00:0.01093)NODE_000000493.00:0.00465)NODE_000000563.00:0.00100;");   starting_tree_12=readTopology("(BHV5:0.20064,(SP1777:0.00316,B589:0.00354)NODE_0000005100.00:0.00660,(216_II:0.00100,(C33:0.00032,(Titanium:0.00000,(C46:0.00000,Cooper:0.00000)NODE_000000299.00:0.00539)NODE_000000167.00:0.00011)NODE_0000000:0.00140)NODE_0000003100.00:0.00416)NODE_0000004100.00:0.00100;");  starting_tree_13=readTopology("(BHV5:0.36810,(SP1777:0.00621,B589:0.00953)NODE_000000576.00:0.00691,(216_II:0.00596,(C33:0.00013,(Titanium:0.00037,(C46:0.00013,Cooper:0.00202)NODE_000000290.00:0.00031)NODE_000000154.00:0.00015)NODE_0000000:0.00485)NODE_0000003100.00:0.00780)NODE_000000497.00:0.00100;"); starting_tree_14=readTopology("(B589:0.00177,BHV5:0.14764,(SP1777:0.00086,(216_II:0.00929,(Cooper:0.00520,(C33:0.00000,(C46:0.00000,Titanium:0.00000)NODE_000000132.00:0.00000)NODE_0000000:0.00029)NODE_000000262.00:0.00732)NODE_0000003100.00:0.00354)NODE_000000483.00:0.00431)NODE_000000583.00:0.00100;");  net_0=readTopology("(((((((BHV5:0.0667816)#H4:0.0540399::0.934345)#H0:0.0838726::0.721054)#H1:1e-06::0.618076,(B589:0.00243892)#H5:0.00247844::0.150925)NODE_000000599.00:0.000967897,(((SP1777:0.00199274,#H5:1e-06::0.849075):0.00142146,#H0:2e-06::0.278946):0.00106455)#H3:1e-06::0.682729):0.000793245,(((216_II:0.000557625,#H4:1e-06::0.0656548):0.0012173,(((((C46:1.13629e-05,(Cooper:0.000504523)#H7:0.000488372::0.309017)NODE_000000287.00:1.38276e-05,#H7:1e-06::0.690983):0.000103885,Titanium:7.99815e-05)NODE_000000172.00:0.000153537,C33:0.000185147)NODE_0000000:0.00010237)#H2:0.000887744::0.731935)NODE_0000003100.00:0.00165082)#H6:1e-06::0.75248):0.00112356,(((#H2:0.00407879::0.268065,#H3:0.00083503::0.317271):3.18362e-05,#H6:0.0018156::0.24752):0.00589199,#H1:0.113802::0.381924):0.00197562)NODE_0000004100.00;");  net_1=readTopology("((((((((B589:0.00235153,(SP1777:0.00196809)#H4:1e-06::0.830798):0.0027904,#H4:0.00150477::0.169202)NODE_000000596.00:1e-06,((((((C46:1.12057e-05,(Cooper:0.000501885)#H5:0.000485768::0.309017)NODE_000000274.00:1.38604e-05,#H5:1e-06::0.690983):0.000103121,Titanium:7.97785e-05)NODE_000000194.00:0.000152737,C33:0.000183986)NODE_0000000:0.000965641,(216_II:1e-06)#H1:0.00165243::0.732512)NODE_000000399.00:0.00261693)#H3:1e-06::0.822789):0.00995,#H1:0.00571558::0.267488):0.0017597,#H3:1e-06::0.177211):0.00797758)#H0:0.0601267::0.701975,(BHV5:0.0165082)#H2:1e-06::0.388659):0.0825265,(#H0:0.0439599::0.298025,#H2:0.0964795::0.611341):0.0738017)NODE_0000004100.00;");  net_2=readTopology("(((BHV5:0.0976493)#H1:0.0659148::0.680401,(((#H1:1e-06::0.319599,(216_II:0.00171691)#H2:1e-06::0.275811):0.00455252,((((((C46:1.12944e-05,(Cooper:0.000503162)#H5:0.000486983::0.309017)NODE_000000283.00:1.38698e-05,#H5:1e-06::0.690983):0.000103331,Titanium:8.01561e-05)NODE_000000169.00:0.000153111,C33:0.000184764)NODE_0000000:0.000956394,#H2:1e-06::0.724189):0.00233281,(((SP1777:0.00192153,(B589:0.00241531)#H4:1e-06::0.793918):0.00195339,#H4:0.0022309::0.206082)NODE_0000005100.00:0.00104678)#H3:1e-06::0.829518):0.00279326)NODE_0000003100.00:1e-06)#H0:0.0593645::0.589172):0.0276427,(#H0:0.00646856::0.410828,#H3:0.0024996::0.170482)NODE_0000004100.00:0.107625)NODE_0000006;");   net_3=readTopology("(((((((((((BHV5:0.03335)#H2:0.128382::0.605345)#H0:0.0482041::0.533902,#H2:1e-06::0.394655):0.0621193)#H1:1e-06::0.749318,(216_II:0.000793618)#H3:0.00113869::0.268355):0.00410144,((((Titanium:7.99148e-05,((C46:1.13626e-05,(Cooper:0.000503603)#H6:0.000488075::0.309017)NODE_0000001100.00:1.39209e-05,#H6:1e-06::0.690983):0.00010402)NODE_000000234.00:0.000153358,C33:0.000185445)NODE_0000000:0.00098999,#H3:0.000875893::0.731645)NODE_0000003100.00:0.00269453)#H4:0.00105073::0.221353):0.00891083,#H4:1e-06::0.778647):0.000545237,(SP1777:0.00200748,(B589:0.00101673)#H5:0.0015024::0.866931)NODE_000000587.00:0.00209412):0.00607325,#H5:1e-06::0.133069):0.00607225,#H0:1e-06::0.466098):0.00607225,#H1:0.0423117::0.250682)NODE_0000004100.00;");   net_4=readTopology("(((((SP1777:0.00183066,(B589:0.0024754)#H3:1e-06::0.729576):0.00166087,#H3:0.00155743::0.270424)NODE_000000547.00:0.000374323,(((((Titanium:8.04985e-05,((C46:1.14136e-05,(Cooper:0.000507005)#H4:0.000490911::0.309017)NODE_000000145.00:1.38176e-05,#H4:1e-06::0.690983):0.000104297)NODE_000000238.00:0.000154312,C33:0.000185665)NODE_0000000:0.0010348,(216_II:0.00164122)#H1:1e-06::0.748398):0.00267111,((BHV5:0.101403)#H0:2e-06::0.329395,#H1:1e-06::0.251602):0.00764238)NODE_0000003100.00:0.00117398)#H2:1e-06::0.779978):0.00426043,#H2:0.00287589::0.220022):0.00101714,#H0:0.181312::0.670605)NODE_0000004100.00;");   net_5=readTopology("(((((((((216_II:0.0016708,(((((C46:1.12754e-05,(Cooper:0.000497365)#H3:0.000499384::0.309017)NODE_000000246.00:1.39047e-05,#H3:1e-06::0.690983):0.000104202,Titanium:8.0068e-05)NODE_000000162.00:0.000153975,C33:0.000185652)NODE_0000000:0.00108957)#H4:1e-06::0.764188):0.0016208)#H1:0.00132731::0.793675,((B589:0.00299296,SP1777:0.00239553)NODE_000000598.00:0.0022661)#H5:1e-06::0.726691):0.00183407,#H4:0.00168183::0.235812)NODE_000000387.00:0.0061756,#H5:0.0016023::0.273309):0.00524404,#H1:1e-06::0.206325):0.0454939)#H0:0.0692778::0.586565,(BHV5:0.0109684)#H2:1e-06::0.194817):0.0539203,(#H2:0.040483::0.805183,#H0:1e-06::0.413435):0.0748925)NODE_0000004100.00;");   net_6=readTopology("(((BHV5:0.0346838)#H2:0.0464512::0.412423,((((((((C46:1.12481e-05,(Cooper:0.000501582)#H5:0.00048511::0.309017)NODE_000000189.00:1.38604e-05,#H5:1e-06::0.690983):0.00010299,Titanium:8.00739e-05)NODE_000000261.00:0.000152559,C33:0.000184659)NODE_0000000:0.000910734,(216_II:0.00180974)#H3:1e-06::0.720701):0.00262375,((SP1777:0.00195605,(B589:0.000873601)#H4:0.00149848::0.864391)NODE_000000580.00:0.00299949)#H1:1e-06::0.830138):0.00544343,#H3:0.00247083::0.279299)NODE_0000003100.00:0.00604451)#H0:1e-06::0.805487):0.000559053,(((#H4:0.00591249::0.135609,#H1:1e-06::0.169862):0.00676287,#H0:1e-06::0.194513):0.079769,#H2:0.00296433::0.587577):0.177922)NODE_000000480.00;");           net_7=readTopology("(((((((216_II:2.1046e-05)#H2:0.00603954::0.268297,((((((C46:1.1271e-05,(Cooper:0.000503824)#H4:0.000487202::0.309017)NODE_000000238.00:1.39445e-05,#H4:1e-06::0.690983):0.000103827,Titanium:7.96343e-05)NODE_000000196.00:0.000153229,C33:0.000184887)NODE_0000000:0.00100749,#H2:0.00162142::0.731703)NODE_000000399.00:0.0030881)#H3:1e-06::0.176004):0.00152845,(((SP1777:0.00198212,(B589:0.0022578)#H5:1e-06::0.899564):0.0022882,#H3:1e-06::0.823996):0.000225479,#H5:0.00798151::0.100436)NODE_000000591.00:0.00790383):0.000678742)#H1:0.0824801::0.777733)#H0:0.0650945::0.582854,#H1:1e-06::0.222267):0.0233901,(BHV5:0.00750113,#H0:1e-06::0.417146):0.135635)NODE_0000004100.00;");   net_8=readTopology("((((((((((Titanium:7.9733e-05,((C46:1.12134e-05,(Cooper:0.000502111)#H5:0.000485789::0.309017)NODE_000000153.00:1.3866e-05,#H5:1e-06::0.690983):0.000103295)NODE_000000243.00:0.000152818,C33:0.000184396)NODE_0000000:0.000985528,(216_II:1e-06)#H0:0.00172429::0.731512)NODE_0000005100.00:0.00268594,((BHV5:1e-06)#H2:0.0613068::0.694298)#H1:0.0956093::0.169398):0.000112274)#H3:1e-06::0.810162,(B589:0.00226526)#H4:0.00345198::0.189372):0.000117368,(SP1777:0.00185582,#H4:1e-06::0.810628):0.0024499)NODE_0000004100.00:0.00921252,(#H0:0.00603239::0.268488,#H3:0.000700342::0.189838):0.000700342)NODE_0000003100.00:0.0896769,#H2:1e-06::0.305702):0.108327,#H1:0.0613068::0.830602)NODE_0000006;");   net_9=readTopology("(((((216_II:8.05607e-05)#H0:2e-06::0.267783,((BHV5:0.0953266)#H2:0.0679223::0.653284)#H1:1e-06::0.287529):0.00183186,#H2:1e-06::0.346716):0.00312621,((((((C46:1.1127e-05,(Cooper:0.000503619)#H5:0.000487603::0.309017)NODE_000000247.00:1.41001e-05,#H5:1e-06::0.690983):0.000103577,Titanium:7.99068e-05)NODE_000000172.00:0.000153317,C33:0.000184945)NODE_0000000:0.000984074,#H0:0.00164646::0.732217)NODE_0000003100.00:0.00274611)#H3:0.00102989::0.21077):0.00249338,(((#H3:1e-06::0.78923,(B589:0.00233288)#H4:0.00549829::0.174058):1e-06,(SP1777:0.00198777,#H4:1e-06::0.825942):0.00254328)NODE_000000596.00:0.00500778,#H1:0.132011::0.712471):0.00183402)NODE_0000004100.00;");   net_10=readTopology("((((BHV5:0.0962818)#H0:1e-06::0.335065,(216_II:0.00167089)#H1:1e-06::0.252698):0.00219009,(#H0:0.109958::0.664935)#H3:1e-06::0.295446):0.00219009,(((((SP1777:0.00184877,(B589:0.00236392)#H4:1e-06::0.67154):0.00142384,#H4:0.00162194::0.32846)NODE_0000004100.00:0.00203377)#H2:0.0014518::0.267013,#H3:0.0935813::0.704554):0.00469577,((((((C46:1.12653e-05,(Cooper:0.000501262)#H5:0.000485194::0.309017)NODE_000000267.00:1.37661e-05,#H5:1e-06::0.690983):0.000102857,Titanium:7.99823e-05)NODE_000000158.00:0.000152513,C33:0.000184518)NODE_0000000:0.000935398,#H1:1e-06::0.747302):0.00207549,#H2:1e-06::0.732987):0.00281932)NODE_0000003100.00:0.00148001)NODE_000000591.00;");   net_11=readTopology("(((((((216_II:0.00171253)#H1:0.00164123::0.272339,(BHV5:0.0929545)#H0:1e-06::0.335733):0.00164123,(((((Titanium:7.98577e-05,((C46:1.12458e-05,(Cooper:0.00050373)#H5:0.000487826::0.309017)NODE_000000268.00:1.39586e-05,#H5:1e-06::0.690983):0.000103665)NODE_000000158.00:0.000153356,C33:0.00018501)NODE_0000000:0.000990777,#H1:1e-06::0.727661):0.00225752,(#H0:0.110159::0.664267)#H3:1e-06::0.197852)NODE_000000493.00:0.000459585)#H2:0.000822853::0.171954)NODE_0000003100.00:0.00752311,#H3:0.0931263::0.802148):0.00288056,#H2:1e-06::0.828046):0.000375615,(SP1777:0.00201091,(B589:0.00230973)#H4:1e-06::0.863954):0.00224299):0.00325433,#H4:0.00271794::0.136046)NODE_000000563.00;");   net_12=readTopology("((C33:0.000184761,(((((216_II:0.000489061,((BHV5:0.0243221)#H0:0.0536742::0.467803)#H1:1e-06::0.137712):0.00131772)#H2:0.00397393::0.285334,((((#H0:0.326028::0.532197,#H1:0.0137441::0.862288):0.0137441,(B589:0.000637203)#H4:0.00441706::0.200336):0.00441606,(SP1777:0.00200709,#H4:0.00150022::0.799664)NODE_0000005100.00:0.00196394)NODE_0000003100.00:0.00277182)#H3:0.0058918::0.27126):0.00447908,#H3:1e-06::0.72874):0.000367031,#H2:1e-06::0.714666):0.000982143)NODE_0000000:5.97825e-05,(((C46:1.12488e-05,(Cooper:0.000500983)#H5:0.000485674::0.309017)NODE_000000299.00:1.38543e-05,#H5:1e-06::0.690983):0.000103404,Titanium:7.95803e-05)NODE_000000167.00:9.25799e-05)NODE_0000004100.00;");  net_13=readTopology("(((((BHV5:0.0807945)#H2:2e-06::0.203297,(216_II:0.00167825)#H0:0.00268754::0.253045):0.00214217,((((((C46:1.16106e-05,(Cooper:0.000516106)#H4:0.000499287::0.309017)NODE_000000290.00:1.40987e-05,#H4:1e-06::0.690983):0.000106188,Titanium:8.18422e-05)NODE_000000154.00:0.000156948,C33:0.000189195)NODE_0000000:0.00102585,#H0:1e-06::0.746955):0.0013827)#H1:0.00191311::0.229363):0.00540156,((SP1777:0.00203967,(B589:0.00240999)#H3:1e-06::0.841226):0.00262726,(#H1:1e-06::0.770637,#H3:0.00371512::0.158774):0.00150306)NODE_000000576.00:0.00422956)NODE_0000003100.00:1e-06,#H2:0.229922::0.796703)NODE_000000497.00;");  net_14=readTopology("((((BHV5:0.0961651)#H0:1e-06::0.33529,(216_II:0.00166521)#H1:1e-06::0.252704):0.00219091,(#H0:0.110431::0.66471)#H3:1e-06::0.294508):0.00219091,(((((SP1777:0.00184998,(B589:0.00236426)#H4:1e-06::0.670798):0.00141692,#H4:0.00161584::0.329202)NODE_000000483.00:0.00204129)#H2:0.00143525::0.267034,#H3:0.0933968::0.705492):0.00469103,((((Titanium:8.00789e-05,((C46:1.12661e-05,(Cooper:0.00050127)#H5:0.000485175::0.309017)NODE_000000262.00:1.3764e-05,#H5:1e-06::0.690983):0.000102853)NODE_000000132.00:0.000152478,C33:0.000184526)NODE_0000000:0.000935954,#H1:1e-06::0.747296):0.00207794,#H2:1e-06::0.732966):0.0027831)NODE_0000003100.00:0.00162558)NODE_000000583.00;")

raw_net_12 = readTopology("((C33:0.000184761,(((((216_II:0.000489061,((BHV5:0.0243221)#H0:0.0536742::0.467803)#H1:1e-06::0.137712):0.00131772)#H2:0.00397393::0.285334,((((#H0:0.326028::0.532197,#H1:0.0137441::0.862288):0.0137441,(B589:0.000637203)#H4:0.00441706::0.200336):0.00441606,(SP1777:0.00200709,#H4:0.00150022::0.799664)NODE_0000005100.00:0.00196394)NODE_0000003100.00:0.00277182)#H3:0.0058918::0.27126):0.00447908,#H3:1e-06::0.72874):0.000367031,#H2:1e-06::0.714666):0.000982143)NODE_0000000:5.97825e-05,(((C46:1.12488e-05,(Cooper:0.000500983)#H5:0.000485674::0.309017)NODE_000000299.00:1.38543e-05,#H5:1e-06::0.690983):0.000103404,Titanium:7.95803e-05)NODE_000000167.00:9.25799e-05)NODE_0000004100.00;"); 

###
# Define vectors of trees/networks
inferred_networks=[net_0, net_1, net_2, net_3, net_4, net_5, net_6, net_7, net_8,net_9, net_10, net_11, net_12, net_13, net_14]

starting_trees=[starting_tree_0, starting_tree_1, starting_tree_2, starting_tree_3, starting_tree_4, starting_tree_5, starting_tree_6, starting_tree_7, starting_tree_8, starting_tree_9, starting_tree_10, starting_tree_11, starting_tree_12, starting_tree_13, starting_tree_14]

# Remove all three cycles from the networks
for network in inferred_networks
    shrink3cycles!(network,false) # why false? I have no idea what that option does.
end

# Assign BHV5 the root for all starting trees
for tree in starting_trees
    rootatnode!(tree,1)
end
# Looking at the plots, we see that trees 10, 11, and 14 are rooted at B589. We
# want all the plots to be rooted at BHV5, so let's do that:
rootatnode!(starting_tree_10,2)
rootatnode!(starting_tree_11,2)
rootatnode!(starting_tree_14,2)

# Cosmetic changes to starting trees
PhyloNetworks.rotate!(starting_tree_10,-2)
PhyloNetworks.rotate!(starting_tree_11,-2)
PhyloNetworks.rotate!(starting_tree_14,-2)

# Cosmetic changes to networks
PhyloNetworks.rootatnode!(net_0,7)
map(i->PhyloNetworks.rotate!(net_0,i),[-3,-22,-24,-15,-11,-12])
map(i->PhyloNetworks.rotate!(net_1,i),[-7,-2,-3,-20])
PhyloNetworks.rotate!(net_2,-10)
PhyloNetworks.rootatnode!(net_3,-8)
map(i->PhyloNetworks.rotate!(net_3,i),[-5,-6,-4,15,-3])
map(i->PhyloNetworks.rotate!(net_4,i),[17,-10])
map(i->PhyloNetworks.rotate!(net_5,i),[-3,-2,-8])
PhyloNetworks.rootatnode!(net_6,-20)
map(i->PhyloNetworks.rotate!(net_6,i),[-3,-7,16,-18,-20])
map(i->PhyloNetworks.rotate!(net_7,i),[-17,11])
PhyloNetworks.rootatnode!(net_8,-3)
map(i->PhyloNetworks.rotate!(net_8,i),[-6,19,-20])
PhyloNetworks.rootatnode!(net_9,19)
map(i->PhyloNetworks.rotate!(net_9,i),[19,-18,-3,-16])
PhyloNetworks.rootatnode!(net_10,-9)
map(i->PhyloNetworks.rotate!(net_10,i),[-14,19,-3,-15])
PhyloNetworks.rootatnode!(net_11,-5)
map(i->PhyloNetworks.rotate!(net_11,i),[-5,-7,-12,14])
PhyloNetworks.rootatnode!(net_12,-13)
PhyloNetworks.rotate!(net_12,10)
PhyloNetworks.rotate!(net_13,18)
map(i->PhyloNetworks.rotate!(net_13,i),[-5,-4,-18])
PhyloNetworks.rootatnode!(net_14,-9)
map(i->PhyloNetworks.rotate!(net_14,i),[19,-15,-3])

# Cosmetic changes to best network with 3-cycles
PhyloNetworks.rootatnode!(raw_net_12,-5)
PhyloNetworks.rotate!(raw_net_12,-5)
PhyloNetworks.rotate!(raw_net_12,-6)
PhyloNetworks.rotate!(raw_net_12,-8)
PhyloNetworks.rotate!(raw_net_12,11)
PhyloNetworks.rotate!(raw_net_12,-13)

# Define vectors of trees/networks (have to run this again)
inferred_networks=[net_0, net_1, net_2, net_3, net_4, net_5, net_6, net_7, net_8,net_9, net_10, net_11, net_12, net_13, net_14]

starting_trees=[starting_tree_0, starting_tree_1, starting_tree_2, starting_tree_3, starting_tree_4, starting_tree_5, starting_tree_6, starting_tree_7, starting_tree_8, starting_tree_9, starting_tree_10, starting_tree_11, starting_tree_12, starting_tree_13, starting_tree_14]
### END CODE BLOCK

#_______________________________________________________________________________
#
# Plot all the output networks
#_______________________________________________________________________________

### START CODE BLOCK
# Setup the layout
R"png(filename='multistart-titanium-all-output-networks.png', width=36, height=24, units='cm', res=300)"
R"layout(matrix(0:15, nrow=4, ncol=4, byrow=TRUE))"
R"layout.show(n=15)" # shows the layout
R"par"(mar=[1,1,1,1]) 

# Plot each network
for k in 0:14
    net=eval(Meta.parse("net_$k")) # There must be better way to do parameter expansion
    PhyloPlots.plot(net, useedgelength=false, showedgenumber=true, shownodenumber=true, showgamma=true)
    R"mtext"("NetRAX Network $k")
end
R"dev.off()"
### END CODE BLOCK

#_______________________________________________________________________________
#
# Plot all the starting trees
#_______________________________________________________________________________

### START CODE BLOCK
R"png(filename='all-start-networks.png', width=36, height=24, units='cm', res=300)"

# Setup the layout
R"layout(matrix(0:15, nrow=4, ncol=4, byrow=TRUE))"
R"layout.show(n=15)" # shows the layout
R"par"(mar=[1,1,1,1]) 

# Plot each tree
for k in 0:14
    tree=eval(Meta.parse("starting_tree_$k")) # There must be better way to do parameter expansion
    PhyloPlots.plot(tree, useedgelength=false, showedgenumber=true, shownodenumber=true, showgamma=true)
    R"mtext"("Starting Tree $k")
end

R"dev.off()"
### END CODE BLOCK

# We observe that there are 6 tree topologies represented by the 15 starting trees:
# Topology 1: 0,1,2,5,7,9,12,13
# Topology 2: 3,4
# Topology 3: 6
# Topology 4: 8
# Topology 5: 10
# Topology 6: 11,14


#_______________________________________________________________________________
#
# Comparing Trees
#_______________________________________________________________________________
### BEGIN CODE BLOCK
map(i->hardwiredClusterDistance(starting_trees[i],inferred_networks[i], false), 1:15)
pcd_output_networks=[]
pcd_starting_trees=[]
for i in 1:15
    for j in (i+1):15
        append!(pcd_output_networks,
                hardwiredClusterDistance(inferred_networks[i],
                                         inferred_networks[j],
                                         false))
        append!(pcd_starting_trees,
                hardwiredClusterDistance(starting_trees[i],
                                         starting_trees[j],
                                         false))
    end
end
x=pcd_starting_trees
y=pcd_output_networks

sum(x)/length(x)
sum(y)/length(y)
length(x)

Pkg.add("PyPlot")
using Plots; pyplot()
Plots.PyPlotBackend()

Pkg.add("Distributions")
using Random, Distributions
bins = 0:25
PyPlot.hist(x, bins, alpha=0.5, label="Starting Trees")
PyPlot.hist(y, bins, alpha=0.5, label="Output Networks")
PyPlot.legend(loc="upper right")
PyPlot.show()
### END CODE BLOCK
y_three_cycles_removed
y_unprocessed


# bins = 0:25
# PyPlot.hist(y_three_cycles_removed, bins, alpha=0.5, label="After Removing 3-cycles")
# PyPlot.hist(y_unprocessed, bins, alpha=0.5, label="Before Removing 3-cycles")
# PyPlot.legend(loc="upper right")
# PyPlot.show()

# histogram(x)
# histogram(y)

# Compare the topology of the starting tree with that of the major tree of the
# corresponding NetRAX inferred network
### BEGIN CODE BLOCK
z=[]
for i in 1:15
    append!(z,hardwiredClusterDistance(majorTree(inferred_networks[i])
,starting_trees[i], false))
end
z
# Results:


### END CODE BLOCK




#_______________________________________________________________________________
#
# Plot the best network and associated tree
# see http://crsl4.github.io/PhyloNetworks.jl/latest/man/snaq_plot/
#_______________________________________________________________________________

# first run the above code to remove 3-cycles, reroot, and to rotate the
# networks

k=13 # the network number to plot. Network 12 was the best  best network for the set1b run. Since netrax starts counting at zero, we take k=13.
tree = starting_trees[k]
net = inferred_networks[k]

# Setup the layout
R"png(filename='set1b_best_network.png', width=36, height=24, units='cm', res=300)"
R"layout(matrix(1, nrow=1, ncol=1, byrow=TRUE))"
R"par"(mar=[1,1,1,1])
PhyloPlots.plot(net, tipcex=1.5, edgewidth=2, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true)
#R"mtext"("Best Netrax network",cex=2)
R"dev.off()"


R"png(filename='set1b_best_network_starting_tree.png', width=36, height=24, units='cm', res=300)"
R"layout(matrix(1, nrow=1, ncol=1, byrow=TRUE))"
R"par"(mar=[1,1,1,1])
PhyloPlots.plot(tree, tipcex=1.5, edgewidth=2,  useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true)
#R"mtext"("Netrax Starting Tree for Best Network Run",cex=2)
R"dev.off()"


R"png(filename='set1b_best_network_with_3_cycles.png', width=36, height=24, units='cm', res=300)"
R"layout(matrix(1, nrow=1, ncol=1, byrow=TRUE))"
R"par"(mar=[1,1,1,1])
PhyloPlots.plot(raw_net_12,  tipcex=1.5, edgewidth=2, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true)
#R"mtext"("Best Netrax network with 3-cycles", cex=2)
R"dev.off()"


