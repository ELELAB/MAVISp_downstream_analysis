python rank_candidates_v2.py \
  -i POLE_filtered_indicators.csv \
  -o POLE_ranked.csv \
  --config rank_config.txt \
  --am-class-col "AlphaMissense classification" \
  --am-class-value pathogenic \
  --demask-col "DeMaSk delta fitness" \
  --demask-threshold 0.25 \
  --mutation-col Mutation \
  --quote-all  
  #--top 50 

