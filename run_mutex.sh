module load python/2.7.6-vanilla

source ./MUTPARAMS.sh



#MAKE COMMAND NAME out of parameters. Then run.
COMMAND="python mutex_triangles.py"

if [ $MUTATIONMATRIX ]
then
    COMMAND=${COMMAND}" -m $MUTATIONMATRIXPATH"
fi
if [ $MINFREQ ]
then
    COMMAND=${COMMAND}" -mf $MINFREQ"
fi

if [ $TOP_PERCENTILE ]
then
COMMAND=${COMMAND}" -tp $TOP_PERCENTILE"
fi

if [ $MINPERCENTILE ]
then
COMMAND=${COMMAND}" -mperc $MINPERCENTILE"
fi

if [ $TOP_NUMBER ]
then
COMMAND=${COMMAND}" -tn $TOP_NUMBER"
fi
if [ $MUTEXPAIRFILE ]
then
COMMAND=${COMMAND}" -mdf $MUTEXPAIRFILE"
fi

if [ $COOCCURPAIRFILE ]
then
COMMAND=${COMMAND}" -cdf $COOCCURPAIRFILE"
fi
if [ $OUTPUTPREFIX ]
then
COMMAND=${COMMAND}" -o $OUTPUTPREFIX"
fi
if [ $PATIENTFILE ]
then
COMMAND=${COMMAND}" -pf $PATIENTFILE"
fi
if [ $GENEFILE ]
then
COMMAND=${COMMAND}" -gf $GENEFILE"
fi
if [ $PATIENTBLACKLIST ]
then
COMMAND=${COMMAND}" -pbf $PATIENTBLACKLIST"
fi
if [ $GENEBLACKLIST ]
then
COMMAND=${COMMAND}" -gbf $GENEBLACKLIST"
fi
if [ $GENELIST1 ]
then
COMMAND=${COMMAND}" -gl1 $GENELIST1"
fi
if [ $GENELIST2 ]
then
COMMAND=${COMMAND}" -gl2 $GENELIST2"
fi

if [ $PAIRLIST ]
then
COMMAND=${COMMAND}" -plf $PAIRLIST"
fi


if [ $NUMPROCESSES ]
then
COMMAND=${COMMAND}" -pcn $NUMPROCESSES"
fi
if [ $GROUPTYPE ]
then
COMMAND=${COMMAND}" -gt $GROUPTYPE"
fi
if [ $PAIRTYPE ]
then
COMMAND=${COMMAND}" -pt $PAIRTYPE"
fi
if [ $TRIPLESTATS ]
then
COMMAND=${COMMAND}" -ctst $TRIPLESTATS"
fi
if [ $TRIPLETSCORES ]
then
COMMAND=${COMMAND}" -ctsc $TRIPLETSCORES"
fi
if [ $LOCALEDGEBETWEENNESS ]
then
COMMAND=${COMMAND}" -leb $LOCALEDGEBETWEENNESS"
fi
if [ $PLOTDISTRIBUTION ]
then
COMMAND=${COMMAND}" -pd $PLOTDISTRIBUTION"
fi
if [ $PLOTPATIENTDISTRIBUTION ]
then
COMMAND=${COMMAND}" -ppd $PLOTPATIENTDISTRIBUTION"
fi

if [ $FILTERCOOCCUR ]
then
COMMAND=${COMMAND}" -fc $FILTERCOOCCUR"
fi

if [ $USEDOWNSTREAM ]
then
COMMAND=${COMMAND}" -ud $USEDOWNSTREAM"
fi

if [ $FILTERCODING ]
then
COMMAND=${COMMAND}" -fcoding $FILTERCODING"
fi

if [ $GAINLOSSCOMBINE ]
then
COMMAND=${COMMAND}" -glc $GAINLOSSCOMBINE"
fi



if [ $FILTERMUTEXGAINLOSS ]
then
COMMAND=${COMMAND}" -fmgl $FILTERMUTEXGAINLOSS"
fi

if [ $FILTERMUTEXSAMESEGMENT ]
then
COMMAND=${COMMAND}" -fmss $FILTERMUTEXSAMESEGMENT"
fi

if [ $FILTERCOOCCURSAMESEGMENT ]
then
COMMAND=${COMMAND}" -fcss $FILTERCOOCCURSAMESEGMENT"
fi

if [ $FCSS_RATIO ]
then
COMMAND=${COMMAND}" -fcss_rt $FCSS_RATIO"
fi

if [ $FCSS_DIFFRATIO ]
then
COMMAND=${COMMAND}" -fcss_mfdrt $FCSS_DIFFRATIO"
fi

if [ $FCSS_COVERAGE ]
then
COMMAND=${COMMAND}" -fcss_ct $FCSS_COVERAGE"
fi

if [ $FCSS_PVALUE ]
then
COMMAND=${COMMAND}" -fcss_pt $FCSS_PVALUE"
fi


if [ $BINGENESBYPAIRS ]
then
COMMAND=${COMMAND}" -bgbp $BINGENESBYPAIRS"
fi


if [ $RATIOPERC ]
then
COMMAND=${COMMAND}" -r_p $RATIOPERC"
fi

if [ $SCOREPERC ]
then
COMMAND=${COMMAND}" -s_p $SCOREPERC"
fi

if [ $COVERAGEPERC ]
then
COMMAND=${COMMAND}" -c_p $COVERAGEPERC"
fi
if [ $COMBSCOREPERC ]
then
COMMAND=${COMMAND}" -cs_p $COMBSCOREPERC"
fi
if [ $MPROB ]
then
COMMAND=${COMMAND}" -mp $MPROB"
fi
if [ $CPROB ]
then
COMMAND=${COMMAND}" -cp $CPROB"
fi
if [ $MAXOVERLAP ]
then
COMMAND=${COMMAND}" -mo $MAXOVERLAP"
fi
if [ $MINCOOCCUR ]
then
COMMAND=${COMMAND}" -mc $MINCOOCCUR"
fi
if [ $MINRATIO ]
then
COMMAND=${COMMAND}" -mcr $MINRATIO"
fi
if [ $MINDISTANCETHRESHOLD ]
then
COMMAND=${COMMAND}" -cdt $MINDISTANCETHRESHOLD"
fi




# Run mutex_triangles.py
echo
echo “**************************************************************************”
echo Running mutex on $PREFIX
echo
echo Command is $COMMAND
$COMMAND

echo
echo $PREFIX finished
echo “**************************************************************************”
echo