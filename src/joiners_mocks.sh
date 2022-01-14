

filtered_mocks=/users/tg/jwillis/SLL/Part_1/SLL1_paper/SRA_submission/mock_communities/filtered_mocks

for mock in HM-782D-plate1  HM-783D-plate1  HM-782D3-plate2 HM-782D4-plate2 HM-783D3-plate2 HM-783D4-plate2 HM-782D-plate3  HM-783D-plate3 \
    HM-782D-plate4  HM-783D-plate4  NTC1-plate4     NTC2-plate4     HM-782D-plate5  HM-783D-plate5  NTC1-plate5     NTC2-plate5 \
    HM-782D-plate6  HM-783D-plate6  NTC1-plate6     NTC2-plate6     HM-782D-plate7  HM-783D-plate7  HM-782D-plate8  HM-783D-plate8 \
    H2Od-plate9     HM-782D-plate9  HM-783D-plate9; do

  echo ""
  echo `date` $mock

  # join by PEAR v0.9.6
  PEAR_dir=$filtered_mocks/joined/PEAR
  pear -f $filtered_mocks/${mock}_F_filt.fastq.gz -r $filtered_mocks/${mock}_R_filt.fastq.gz -o $PEAR_dir/${mock}_PEAR
  gzip -f $PEAR_dir/${mock}_PEAR.assembled.fastq
  rm $PEAR_dir/${mock}_PEAR.discarded.fastq
  rm $PEAR_dir/${mock}_PEAR.unassembled.*.fastq

  echo ""

  # join by pandaseq v2.11 with -A simple_bayesian (option for algorithm)
  pandaseq_sb=$filtered_mocks/joined/pandaseq_sb
  pandaseq -f $filtered_mocks/${mock}_F_filt.fastq.gz -r $filtered_mocks/${mock}_R_filt.fastq.gz -F -A simple_bayesian \
    -w $pandaseq_sb/${mock}_pandaseq.simple_bayesian.fastq 2> $pandaseq_sb/${mock}_pandaseq.simple_bayesian.fastq.log
  gzip -f $pandaseq_sb/${mock}_pandaseq.simple_bayesian.fastq

  echo ""

  # join by pandaseq v2.11 with -A ea_util (option for algorithm) (fastqJoin)
  pandaseq_ea=$filtered_mocks/joined/pandaseq_ea
  pandaseq -f $filtered_mocks/${mock}_F_filt.fastq.gz -r $filtered_mocks/${mock}_R_filt.fastq.gz -F -A ea_util \
    -w $pandaseq_ea/${mock}_pandaseq.ea_util.fastq 2> $pandaseq_ea/${mock}_pandaseq.ea_util.fastq.log
  gzip -f $pandaseq_ea/${mock}_pandaseq.ea_util.fastq

  echo ""

  # join by pandaseq v2.11 with -A flash (option for algorithm)
  pandaseq_flash=$filtered_mocks/joined/pandaseq_flash
  pandaseq -f $filtered_mocks/${mock}_F_filt.fastq.gz -r $filtered_mocks/${mock}_R_filt.fastq.gz -F -A flash \
    -w $pandaseq_flash/${mock}_pandaseq.flash.fastq 2> $pandaseq_flash/${mock}_pandaseq.flash.fastq.log
  gzip -f $pandaseq_flash/${mock}_pandaseq.flash.fastq

  echo ""

  # join by pandaseq v2.11 with -A pear (option for algorithm)
  pandaseq_pear=$filtered_mocks/joined/pandaseq_pear
  pandaseq -f $filtered_mocks/${mock}_F_filt.fastq.gz -r $filtered_mocks/${mock}_R_filt.fastq.gz -F -A pear \
    -w $pandaseq_pear/${mock}_pandaseq.pear.fastq 2> $pandaseq_pear/${mock}_pandaseq.pear.fastq.log
  gzip -f $pandaseq_pear/${mock}_pandaseq.pear.fastq

  echo ""

  # join by pandaseq v2.11 with -A rdp_mle (option for algorithm)
  pandaseq_rdp_mle=$filtered_mocks/joined/pandaseq_rdp_mle
  pandaseq -f $filtered_mocks/${mock}_F_filt.fastq.gz -r $filtered_mocks/${mock}_R_filt.fastq.gz -F -A rdp_mle \
    -w $pandaseq_rdp_mle/${mock}_pandaseq.rdp_mle.fastq 2> $pandaseq_rdp_mle/${mock}_pandaseq.rdp_mle.fastq.log
  gzip -f $pandaseq_rdp_mle/${mock}_pandaseq.rdp_mle.fastq

  echo ""

  # join by pandaseq v2.11 with -A stitch (option for algorithm)
  pandaseq_stitch=$filtered_mocks/joined/pandaseq_stitch
  pandaseq -f $filtered_mocks/${mock}_F_filt.fastq.gz -r $filtered_mocks/${mock}_R_filt.fastq.gz -F -A stitch \
    -w $pandaseq_stitch/${mock}_pandaseq.stitch.fastq 2> $pandaseq_stitch/${mock}_pandaseq.stitch.fastq.log
  gzip -f $pandaseq_stitch/${mock}_pandaseq.stitch.fastq

  echo ""

  # join by pandaseq v2.11 with -A uparse (option for algorithm)
  pandaseq_uparse=$filtered_mocks/joined/pandaseq_uparse
  pandaseq -f $filtered_mocks/${mock}_F_filt.fastq.gz -r $filtered_mocks/${mock}_R_filt.fastq.gz -F -A uparse \
    -w $pandaseq_uparse/${mock}_pandaseq.uparse.fastq 2> $pandaseq_uparse/${mock}_pandaseq.uparse.fastq.log
  gzip -f $pandaseq_uparse/${mock}_pandaseq.uparse.fastq

  echo ""

  # join by fastq-join v1.3.1
  fqj_dir=$filtered_mocks/joined/fqj
  fastq-join $filtered_mocks/${mock}_F_filt.fastq.gz $filtered_mocks/${mock}_R_filt.fastq.gz -o $fqj_dir/${mock}_fqj.
  mv $fqj_dir/${mock}_fqj.join $fqj_dir/${mock}_fqj.join.fastq
  rm $fqj_dir/${mock}_fqj.un*
  gzip -f $fqj_dir/${mock}_fqj.join.fastq

  echo ""

done




