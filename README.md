# SaCas9_off-target_activity
with https://github.com/editasmedicine/pairwise-library-screen

1. find off-target candidate by Cas-Offinder
    result form : 
        WT(guide)  chromosome  location    mis_match_seq   strand  #_of_mis_match

2-1. filter out if #_mis_match >= 2  ==> list of filterd seq by guide(21 bp)
2-2. filter out if #_mis_match == 0 doesn't exist ==> list of seq that #_mis_match == 0 doesn't exist
    result :
        WT(guide)[:21]   #_of_mis_match  #_seq  tot_seq
        AGCT...CGT  0   2   15
        AGCT...AGT  0   1   8

2. get off-target score by predict_activity_single
    input : WT[:21] mis_match_seq[:21]
    
    

####################################################################################    
1) guide RNAs design
    1-1) Initial input file: 201130_CCDS_human_current.txt, 201130_CCDS_mouse_current.txt
    1-2) Input filter
        1-2-1) ccds_id: 길이가 가장 짧은 cds에 해당하는 ccds_id만 사용
        1-2-2) ccds_status: Public만 사용
        1-2-3) match_type: Identical만 사용
    1-3) Reference genome (2200): /media/backup/ref/Ensemble_GRCh38_p13, /media/backup/ref/Ensemble_GRCm38_p6
    1-4) guide RNA 기준
        1-4-1) cds의 5-65% 구간에 cleavage site가 위치하는 guide RNA context를 genomic region에서 찾음 (cleavage site가 junction에 겹치면 제거)
        1-4-2) Cas9별로 길이 및 PAM 조건은 기존에 전달해드린 것과 같아 따로 적지 않았으나 필요하시면 이 파일에 업데이트 하겠습니다.
    1-5)  Output file
        1-5-1) Filtered input
        1-5-2) sgRNA with context sequence (location, strand, cds 내에서 cleavage site의 위치 정보 등을 포함)

2) Cas-offinder (제가 돌릴 부분이라 생략하겠습니다.)

3) SaCas9 off-target activity
    3-1) Initial input file directory: /extdata1/JaeWoo/Project/YoungGwang_CasOffFinder/Output/FirstResult/ 
    3-2) Input filter: mismatch가 0인 target의 개수가 1인 guide RNA만 사용
    3-3) Soff-target:  https://github.com/editasmedicine/pairwise-library-screen 의 predict_activity_single.py 
    3-4) Sguide = 100/(1+ ∑_(𝑖=1)^𝑛▒〖𝑆𝑜𝑓𝑓−𝑡𝑎𝑟𝑔𝑒𝑡(𝑖)〗) 
                = 100/(1 + sum(Soff_pred_activity))
    3-5) Output file
        3-5-1) 각각의 guide별로 mismatch 개수에 따른 target 개수 정리 (e.g. guide RNA sequence, # of 0-bp mismatched targets, ...,  # of 3-bp mismatched targets)
        3-5-2) 각각의 guide 별로 3-4)에서 구한 최종 score
        3-5-3) 3-5-1)과 3-5-2)는 함께 주시면 정리하기 좋겠지만 파일이 너무 커지는 등의 문제가 발생한다면 따로 주셔도 무방합니다.
