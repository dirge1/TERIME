function [I,V]=IVload
% I=[0.00434810905045534,0.00434096133844389,0.00433379747037609,0.00432661636898789,0.00431941688559253,0.00431219779540136,0.00430495779254607,0.00429769548478396,0.00429040938786759,0.00428309791955879,0.00427575939326611,0.00426839201128399,0.00426099385761051,0.00425356289031997,0.00424609693346531,0.00423859366848419,0.00423105062508201,0.00422346517156364,0.00421583450458503,0.00420815563829475,0.00420042539283500,0.00419264038217046,0.00418479700121309,0.00417689141221042,0.00416891953036450,0.00416087700864858,0.00415275922178895,0.00414456124937949,0.00413627785809754,0.00412790348299068,0.00411943220780579,0.00411085774433372,0.00410217341074587,0.00409337210890206,0.00408444630061359,0.00407538798285013,0.00406618866188507,0.00405683932638116,0.00404733041942616,0.00403765180953788,0.00402779276066872,0.00401774190125231,0.00400748719234869,0.00399701589496034,0.00398631453660898,0.00397536887728262,0.00396416387488388,0.00395268365033444,0.00394091145251619,0.00392882962325773,0.00391641956260462,0.00390366169464410,0.00389053543418832,0.00387701915465566,0.00386309015752589,0.00384872464378231,0.00383389768779086,0.00381858321410360,0.00380275397770879,0.00378638154828350,0.00376943629903449,0.00375188740073824,0.00373370282161109,0.00371484933365217,0.00369529252610582,0.00367499682668316,0.00365392553116386,0.00363204084196731,0.00360930391623533,0.00358567492390643,0.00356111311618166,0.00353557690468548,0.00350902395151096,0.00348141127020673,0.00345269533761626,0.00342283221631776,0.00339177768723931,0.00335948739184126,0.00332591698306948,0.00329102228409461,0.00325475945366738,0.00321708515674466,0.00317795673888002,0.00313733240273217,0.00309517138492928,0.00305143413144259,0.00300608246957163,0.00295907977463024,0.00291039112944845,0.00285998347487104,0.00280782574953912,0.00275388901738338,0.00269814658143464,0.00264057408276302,0.00258114958358720,0.00251985363384287,0.00245666932075730,0.00239158230123911,0.00232458081715054,0.00225565569377786,0.00218480032204804,0.00211201062525068,0.00203728501121010,0.00196062431100917,0.00188203170549341,0.00180151264087855,0.00171907473484854,0.00163472767456457,0.00154848310801078,0.00146035453008146,0.00137035716477140,0.00127850784476789,0.00118482488966420,0.00108932798392271,0.000992038055615863,0.000892977156866453,0.000792168346800004,0.000689635577712059,0.000585403585046032,0.000479497781673167,0.000371944156867858,0.000262769180279345,0.000151999711116078,3.96629126818652e-05];
% V=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000,0.710000000000000,0.720000000000000,0.730000000000000,0.740000000000000,0.750000000000000,0.760000000000000,0.770000000000000,0.780000000000000,0.790000000000000,0.800000000000000,0.810000000000000,0.820000000000000,0.830000000000000,0.840000000000000,0.850000000000000,0.860000000000000,0.870000000000000,0.880000000000000,0.890000000000000,0.900000000000000,0.910000000000000,0.920000000000000,0.930000000000000,0.940000000000000,0.950000000000000,0.960000000000000,0.970000000000000,0.980000000000000,0.990000000000000,1,1.01000000000000,1.02000000000000,1.03000000000000,1.04000000000000,1.05000000000000,1.06000000000000,1.07000000000000,1.08000000000000,1.09000000000000,1.10000000000000,1.11000000000000,1.12000000000000,1.13000000000000,1.14000000000000,1.15000000000000,1.16000000000000,1.17000000000000,1.18000000000000,1.19000000000000,1.20000000000000,1.21000000000000,1.22000000000000,1.23000000000000];
V=[-0.2057
-0.1291
-0.0588
0.0057
0.0646
0.1185
0.1678
0.2132
0.2545
0.2924
0.3269
0.3585
0.3873
0.4137
0.4373
0.459
0.4784
0.496
0.5119
0.5265
0.5398
0.5521
0.5633
0.5736
0.5833
0.59
]';
I=[0.764
0.762
0.7605
0.7605
0.76
0.759
0.757
0.757
0.7555
0.754
0.7505
0.7465
0.7385
0.728
0.7065
0.6755
0.632
0.573
0.499
0.413
0.3165
0.212
0.1035
-0.0100
-0.1230
-0.2100
]';
end
