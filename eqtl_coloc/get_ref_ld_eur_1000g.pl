#! /usr/bin/perl

# Based on https://github.com/cotsapaslab/jlim/blob/master/bin/fetch.refld0.EUR.pl
# Modified to include all unrelated EUR 1000G indiviudals

use strict;
use warnings;

my $vcfdir = $ARGV[0];
my $idxSNPfile = $ARGV[1];
my $outdir = $ARGV[2];

unless ( -d "$outdir" ) {
    `mkdir $outdir`;
}

open(IDX, "<$idxSNPfile") or die "Cannot open $idxSNPfile";

foreach my $line (<IDX>) {

    # header line
    next if $line =~ /^CHR/;

    chomp $line;
    my @window_info = split /\t/, $line;

    my $chrom = $window_info[0];
    my $from = $window_info[3];
    my $to = $window_info[4];

    if ( -f "$outdir/locus.$chrom.$from.$to.txt.gz" ) {
	print "Already exist: $outdir/locus.$chrom.$from.$to.txt.gz\n";
	next;
    }

    print "'$chrom' '$from' '$to'\n";

    `vcf-query $vcfdir/ALL.chr$chrom.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz $chrom:$from-$to -c HG00096,HG00097,HG00099,HG00100,HG00101,HG00102,HG00103,HG00104,HG00105,HG00106,HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116,HG00117,HG00118,HG00119,HG00120,HG00121,HG00122,HG00123,HG00125,HG00126,HG00127,HG00128,HG00129,HG00130,HG00131,HG00132,HG00133,HG00134,HG00135,HG00136,HG00137,HG00138,HG00139,HG00140,HG00141,HG00142,HG00143,HG00145,HG00146,HG00148,HG00149,HG00150,HG00151,HG00152,HG00154,HG00155,HG00156,HG00157,HG00158,HG00159,HG00160,HG00171,HG00173,HG00174,HG00176,HG00177,HG00178,HG00179,HG00180,HG00181,HG00182,HG00183,HG00185,HG00186,HG00187,HG00188,HG00189,HG00190,HG00231,HG00232,HG00233,HG00234,HG00235,HG00236,HG00237,HG00238,HG00239,HG00240,HG00242,HG00243,HG00244,HG00245,HG00246,HG00249,HG00250,HG00251,HG00252,HG00253,HG00254,HG00255,HG00256,HG00257,HG00258,HG00259,HG00260,HG00261,HG00262,HG00263,HG00264,HG00265,HG00266,HG00267,HG00268,HG00269,HG00270,HG00271,HG00272,HG00273,HG00274,HG00275,HG00276,HG00277,HG00278,HG00280,HG00281,HG00282,HG00284,HG00285,HG00288,HG00290,HG00302,HG00303,HG00304,HG00306,HG00308,HG00309,HG00310,HG00311,HG00312,HG00313,HG00315,HG00318,HG00319,HG00320,HG00321,HG00323,HG00324,HG00325,HG00326,HG00327,HG00328,HG00329,HG00330,HG00331,HG00332,HG00334,HG00335,HG00336,HG00337,HG00338,HG00339,HG00341,HG00342,HG00343,HG00344,HG00345,HG00346,HG00349,HG00350,HG00351,HG00353,HG00355,HG00356,HG00357,HG00358,HG00359,HG00360,HG00361,HG00362,HG00364,HG00365,HG00366,HG00367,HG00368,HG00369,HG00371,HG00372,HG00373,HG00375,HG00376,HG00377,HG00378,HG00379,HG00380,HG00381,HG00382,HG00383,HG00384,HG01334,HG01500,HG01501,HG01503,HG01504,HG01506,HG01507,HG01509,HG01510,HG01512,HG01513,HG01515,HG01516,HG01518,HG01519,HG01521,HG01522,HG01524,HG01525,HG01527,HG01528,HG01530,HG01531,HG01536,HG01537,HG01602,HG01603,HG01605,HG01606,HG01607,HG01608,HG01610,HG01612,HG01613,HG01615,HG01617,HG01618,HG01619,HG01620,HG01623,HG01624,HG01625,HG01626,HG01628,HG01630,HG01631,HG01632,HG01668,HG01669,HG01670,HG01672,HG01673,HG01675,HG01676,HG01678,HG01679,HG01680,HG01682,HG01684,HG01685,HG01686,HG01694,HG01695,HG01697,HG01699,HG01700,HG01702,HG01704,HG01705,HG01707,HG01708,HG01709,HG01710,HG01746,HG01747,HG01756,HG01757,HG01761,HG01762,HG01765,HG01766,HG01767,HG01768,HG01770,HG01771,HG01773,HG01775,HG01776,HG01777,HG01779,HG01781,HG01783,HG01784,HG01785,HG01786,HG01789,HG01790,HG01791,HG02215,HG02219,HG02220,HG02221,HG02223,HG02224,HG02230,HG02231,HG02232,HG02233,HG02235,HG02236,HG02238,HG02239,HG04301,HG04302,HG04303,NA06984,NA06985,NA06986,NA06989,NA06994,NA07000,NA07037,NA07048,NA07051,NA07056,NA07347,NA07357,NA10847,NA10851,NA11829,NA11830,NA11831,NA11832,NA11840,NA11843,NA11881,NA11892,NA11893,NA11894,NA11918,NA11919,NA11920,NA11930,NA11931,NA11932,NA11933,NA11992,NA11994,NA11995,NA12003,NA12004,NA12005,NA12006,NA12043,NA12044,NA12045,NA12046,NA12058,NA12144,NA12154,NA12155,NA12156,NA12234,NA12249,NA12272,NA12273,NA12275,NA12282,NA12283,NA12286,NA12287,NA12340,NA12341,NA12342,NA12347,NA12348,NA12383,NA12399,NA12400,NA12413,NA12414,NA12489,NA12546,NA12716,NA12717,NA12718,NA12748,NA12749,NA12750,NA12751,NA12760,NA12761,NA12762,NA12763,NA12775,NA12776,NA12777,NA12778,NA12812,NA12813,NA12814,NA12815,NA12827,NA12828,NA12829,NA12830,NA12842,NA12843,NA12872,NA12873,NA12874,NA12878,NA12889,NA12890,NA20502,NA20503,NA20504,NA20505,NA20506,NA20507,NA20508,NA20509,NA20510,NA20511,NA20512,NA20513,NA20514,NA20515,NA20516,NA20517,NA20518,NA20519,NA20520,NA20521,NA20522,NA20524,NA20525,NA20527,NA20528,NA20529,NA20530,NA20531,NA20532,NA20533,NA20534,NA20535,NA20536,NA20537,NA20538,NA20539,NA20540,NA20541,NA20542,NA20543,NA20544,NA20581,NA20582,NA20585,NA20586,NA20587,NA20588,NA20589,NA20752,NA20753,NA20754,NA20755,NA20756,NA20757,NA20758,NA20759,NA20760,NA20761,NA20762,NA20763,NA20764,NA20765,NA20766,NA20767,NA20768,NA20769,NA20770,NA20771,NA20772,NA20773,NA20774,NA20775,NA20778,NA20783,NA20785,NA20786,NA20787,NA20790,NA20792,NA20795,NA20796,NA20797,NA20798,NA20799,NA20800,NA20801,NA20802,NA20803,NA20804,NA20805,NA20806,NA20807,NA20808,NA20809,NA20810,NA20811,NA20812,NA20813,NA20814,NA20815,NA20816,NA20818,NA20819,NA20821,NA20822,NA20826,NA20827,NA20828,NA20829,NA20831,NA20832 -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' | sed 's/|/\t/g' | gzip > $outdir/locus.$chrom.$from.$to.txt.gz`;
    if ($? != 0) {
	# fail
	unlink "$outdir/locus.$chrom.$from.$to.txt.gz";
	die "failed";
    }

}

close(IDX);

print "Done\n";
exit 0;
