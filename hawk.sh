#!/bin/bash
khufu_dir="/cluster/projects/khufu/qtl_seq_II/khufu_II"
source "$khufu_dir"/utilities/load_modules.sh
##################
helpFunc()
{
   "$khufu_dir"/utilities/logo.sh
   echo -e "
   Usage:              \033[46m$0 -bed asd.bed -panmap asd.panmap -background "asd1,asd2"\033[0m
   \033[46m-/--t\033[0m               number of threads
   \033[46m-/--o\033[0m               output file name
   \033[46m-/--bed\033[0m             bed file: chr\tpos1\tpos2\tparent,trait
   \033[46m-/--panmap\033[0m          input panmap file
   \033[46m-/--background\033[0m      string of background parents in comma delimted format, default no background parent, All to use all parents of the panmap file
   \033[46m-/--MinDep\033[0m          int of the minimum number of variants in region/sample without missing,to be calcuated, default is 5
   \033[46m-/--dominance\033[0m       float (NA or 0:1); NA means het will not be counted, default is partial domimance 1
   \033[46m-/--clean\033[0m           clean tmp folder/files
   \033[46m-/--graph\033[0m  if the run will produce RDs object for visualization [0 or 1; default=0]
   "
   exit 1
}
####################################################################################
t=4
MinDep=5
dominance=1
clean=1
graph=0
####################################################################################
SOPT='t:o:h'
LOPT=('bed' 'panmap' 'graph' 'background' 'MinDep' 'dominance' 'clean')
OPTS=$(getopt -q -a --options ${SOPT} --longoptions "$(printf "%s:," "${LOPT[@]}")" --name "$(basename "$0")" -- "$@")
eval set -- $OPTS
while [[ $# > 0 ]]; do
    case ${1} in
      -h) helpFunc ;;
      -t) t=$2 && shift ;;
      -o) OutFile=$2 && shift ;;
      --bed) bed=$2 && shift ;;
      --panmap) panmap=$2 && shift ;;
      --background) BACK=$2 && shift ;;
      --dominance) dominance=$2 && shift ;;
      --MinDep) MinDep=$2 && shift ;;
      --graph) graph=$2 && shift ;;
      --clean) clean=$2 && shift ;;
      esac
    shift
done
####################################################################################
echo "t=$t"
if [[ -z "$OutFile" ]]; then OutFile=$(echo $panmap | sed "s:.*/::g" | sed "s:$:.hawk:g"); fi
if [[ $bed == "" ]]; then echo "bed file should be provided" ; exit 0; fi
if [[ $panmap == "" ]]; then echo "panmap file should be provided" ; exit 0; fi
echo "bed=$bed"
echo "panmap=$panmap"
echo "BACK=$BACK"
echo "MinDep=$MinDep"
DomCheck=$(echo $dominance | awk '{if($0 == "NA") {print "het loci will not be counted"} else if($0 >= 0 && $0 <= 1){print "dominance = "$0 }else{print "ERROR" } }')
if [[ $DomCheck == "ERROR" ]]; then echo "ERROR: dominance should be betweeen 0 and 1"; exit 0 ; else echo $DomCheck; fi
if [[ $graph == 0 ]]; then echo "# Graph production is disabled"; elif [[ $graph == 1 ]]; then echo "# Graph will be produced, run may take a longer time"; else echo "graph should be 1 or 0"; exit 0; fi
if [[ $clean == 1 ]]; then echo "# tmp folder/files will be removed" ; fi
echo "output=$OutFile"
####################################################################################
#bed="QTLs.bed"
#panmap="KHU124_min2_Smiss1_miss1_maf0.0.panmap"
#dominance=1
#MinDep=5
#BACK="TifNV,CC477"
#t=4
####################################################################################
tmpDir=$(mktemp -d "KhufuEnviron.XXXXXXXXX")
###
cat $bed | tr ',' '\t' > "$tmpDir"/bed.txt
for y in $(seq 1 $(cat $bed | wc -l))
do
   (
   QTLs=($(cat "$tmpDir"/bed.txt | sed -n "$y p" | tr '\t' ' ' ))
   echo "$y"" out of " $(cat "$tmpDir"/bed.txt | wc -l) "${QTLs[@]}"
   chr=${QTLs[0]}
   pos1=${QTLs[1]}
   pos2=${QTLs[2]}
   par=${QTLs[3]}
   qtl=${QTLs[4]}
   cat $panmap |  awk -v chr=$chr -v pos1=$pos1 -v pos2=$pos2 '{if(NR == 1 || ($1==chr && $2 >= pos1 && $2 <= pos2) ) print $0}' > "$tmpDir"/panmap.$y.txt

inx=$(cat $panmap | cut -f 4  | head -1 | tr ',' '\n' | grep -nE $(echo $par | sed "s:;:|:g" ) | cut -d':' -f 1)
# drop not unique or have 0
cat "$tmpDir"/panmap.$y.txt | cut -f 4  | sed "1d" | awk -v inx=$(echo $inx | tr ' ' ',') '{nA=split($0,A,",");nB=split(inx,B,","); X1=""; for(b=1;b<=nB;++b){X1=X1","A[B[b]]}; X1=substr(X1,2); print X1"\t"$0 }' | awk '{if($1 ~ "0") {print $0"\tOUT"} else {print $0"\tIN"} } ' | awk '{nA=split($1,A,","); for (a=1;a<=nA;++a){gsub(A[a],"U",$0) } }'1  | awk '{nA=split($1,A,","); print gsub("U","",$2)"\t"nA"\t"$3}' | awk '{if($1==$2 && $3 == "IN") {print 1} else {print 0} }' | sed "1i1\t1"| paste - "$tmpDir"/panmap.$y.txt  | awk '{if($1==1) print $0}' | cut -f 2- > "$tmpDir"/panmap.$y.txt2A
# drop multiple allele for the parents
cat "$tmpDir"/panmap.$y.txt2A | cut -f 4  | sed "1d" | awk -v inx=$(echo $inx | tr ' ' ',') '{nA=split($0,A,",");nB=split(inx,B,","); X1=""; for(b=1;b<=nB;++b){X1=X1","A[B[b]]}; X1=substr(X1,2); print X1 }' | SortUniqueValuesDS | sed "1i1\t1"| paste - "$tmpDir"/panmap.$y.txt2A  | awk '{nA=split($1,A,","); if(nA==1) print $0}' | sed "1 s:^.*chr:1\tchr:g" > "$tmpDir"/panmap.$y.txt2

   if [[ "$graph" == 1 ]]
   then
   
   cat "$tmpDir"/panmap.$y.txt2 | cut -f 2- | awk '{if(NR==1) {split($4,Par,",");split($0,IDs,"\t")} else{split($3,LEN,",");nP=split($4,P,","); for (i=5;i<=NF;++i){ if($i != "-"){ nA=split($i,A,","); for(a=1;a<=nA;++a) { for(p=1;p<=nP;++p){ if(P[p]==A[a]){print IDs[i],$1,$2,Par[p],LEN[A[a]] }  }} }  }  } }' | awk -v par=$par '{nA=split(par,A,";"); X=0; for(a=1;a<=nA;++a){ if(A[a]==$4){X=1} }; print $1,$2,$3,X,$5  }'  | tr ' ' '\t' | sort -k4,4n | sed "1i id\tchr\tpos\tcall\tlen" | Rscript -e '{args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=TRUE); A$call = as.factor(A$call); library(ggplot2); g1 = ggplot(A,aes(pos,id)) + geom_point(aes(color=call,size=len,alpha=0.5)) + scale_color_manual(values=c("light blue","dark blue")) + ggtitle(gsub("_"," ",args[2])); pdf(args[1],width=42,height=round(length(unique(A$id))/2)); print(g1); dev.off() }' "$tmpDir"/"$OutFile"_"$qtl"_"$par".pdf $(echo ${QTLs[@]} | tr ' ' '_' )

   fi
      ###
   tcount=1
   for x in $(seq 6 $(cat "$tmpDir"/panmap.$y.txt2 | head -1 | wc -w))
   do
      (
      echo "$y.$qtl:  $x out of $(cat "$tmpDir"/panmap.$y.txt2 | head -1 | wc -w)"
      id=$(cat "$tmpDir"/panmap.$y.txt2 | head -1 | tr '\t' '\n' | sed -n "$x p" )
      cat "$tmpDir"/panmap.$y.txt2 | cut -f 1-3,"$x" | sed "1d" | grep -v "-"  > "$tmpDir"/panmap.$y.ind_$x.txt2
      cat "$tmpDir"/panmap.$y.ind_$x.txt2  | grep -v "-" | awk -v dominance=$dominance '{isdominancematch=0; nA=split($4,A,",");  if(nA==1){ if($1==$4){X=1} else {X=0} } else {X=0; for(a=1;a<=nA;++a) { if(A[a]==$1){X=dominance; isdominancematch=1; break}  } }; print $1"\t"$2"\t"$3"\t"X"\t"isdominancematch  }' | grep -vw "NA$" | sed "s:$:\t$id:g" > "$tmpDir"/panmap.$y.ind_$x.txt3
      cat "$tmpDir"/panmap.$y.ind_$x.txt3  | awk 'BEGIN{X1=0; X2=0; ; Sim=0;Cnt=0} {if(NR==1){X1=$3}; if(NR==NR){X2=$3}; Sim=Sim+$4; Cnt+=1  }END{if(Cnt==0){sim=0}else{sim=Sim/Cnt} ;printf "%0.2f\t%s\t%s\t%s\n", sim,Cnt,X1,X2}' | sed "s:^:$x\t$y\t$id\t$qtl\t:g" > "$tmpDir"/"$x"."$y".txt4
      ) &
      if [[ "$tcount" == $t ]]; then  wait; tcount=1; else tcount=$((tcount+1)); fi
   done
   wait
   if [[ "$graph" == 1 ]]
   then
      cat "$tmpDir"/panmap.$y.ind_*.txt3 | awk '{printf "%s\t%s\t%s\t%s\n", $6,$2,$3,$4}' | sed "s:$:\t$qtl:g" > "$tmpDir"/Fore."$y".graph
   fi
   )
done
####################################################################################
cat $tmpDir/*.txt4 | sort -k1,1n -k2,2n  | cut -f 3- > "$tmpDir"/asd.txt4
cat "$tmpDir"/asd.txt4 | awk -v MinDep=$MinDep '{if($4<MinDep){$3="NA"} }'1 > "$tmpDir"/asd.txt5
cat "$tmpDir"/asd.txt5 | awk '{ if(id==$1){ len=len","$6-$5; range=range","$5"-"$6; if($3 != "NA"){S+=$3; N+=1};Sim=Sim","$3; numVar=numVar","$4  } else {Aver=0; if(N>0){Aver=S/N}; printf "%s\t%0.2f\t%s\t%s\t%s\t%s\n", id,Aver,Sim,numVar,range,len  ;id=$1; len=$6-$5; range=$5"-"$6; S=$3; if($3 == "NA"){N=0}else(N=1); Sim=$3; numVar=$4}  }END{ printf "%s\t%0.2f\t%s\t%s\t%s\t%s\n", id,Aver,Sim,numVar,range,len }' | sed "1d" | sort -k2,2nr | sed "1iID\tAverageSimilarity\tSimilarity\tNumberOfVariants\tRange\tlen" > "$tmpDir"/asd.txt7
if [[ "$graph" == 1 ]];then cat "$tmpDir"/Fore.*.graph | sed "1iID\tchr\tpos\tallele\tqtl" > "$tmpDir"/ForeGraph.graph; pdfunite "$tmpDir"/"$OutFile"_*_*.pdf "$OutFile"_ForeHawk.pdf  ;fi
####################################################################################
## background
if [[ -z "$BACK" ]]
then
   cp "$tmpDir"/asd.txt7 $OutFile
   if [[ "$graph" == 1 ]];then cp "$tmpDir"/ForeGraph.graph "$tmpDir"/FullGraph.graph; fi
else
   if [[ $BACK == "All" ]]; then BACK=$(cat $panmap | head -1 | cut -f 4); fi
   ##
   cat "$tmpDir"/panmap.*.txt  | cut -f 1-2 | sort -k1,1 -k2,2n | tr '\t' '_' | grep -v "chr_pos" | sed "s:$:\tF:g" >  "$tmpDir"/Fore.ids
   cat $panmap | cut -f 1-2 | tr '\t' '_'  > "$tmpDir"/All.ids
   merge "$tmpDir"/All.ids  "$tmpDir"/Fore.ids | awk '{if($2=="NA") {print 1} else {print 0} }' | paste - $panmap | awk '{if($1==1) print $0}' | cut -f 2- >  "$tmpDir"/back.panmap
   for back in $(echo $BACK | tr ',' '\n' )
   do
      (
      inx=$(cat $panmap | cut -f 4  | head -1 | tr ',' '\n' | grep -n $back | cut -d':' -f 1)
      cat "$tmpDir"/back.panmap | cut -f 4  | sed "1d" | awk -v inx=$inx '{nA=split($0,A,","); X=0; if(A[inx] == 0){X=0} else { for(a=1;a<=nA;++a){if(A[inx]==A[a]){X+=1} }  }; print X"\t"A[inx] }' | sed "1i1\t1"| paste - "$tmpDir"/back.panmap | awk '{if($1==1) print $0}' | cut -f 2- > "$tmpDir"/back_"$back".txt
      tcount=1
      for x in $(seq 6 $(cat "$tmpDir"/back_"$back".txt | head -1 | wc -w) )
      do
         (
         id=$(cat "$tmpDir"/back_"$back".txt | head -1 | tr '\t' '\n' | sed -n "$x p" )
         echo "BG:  $x out of $(cat "$tmpDir"/back_"$back".txt | head -1 | wc -w)"
         cat "$tmpDir"/back_"$back".txt | cut -f 1-3,$x | grep -v "-" | awk -v dominance=$dominance '{isdominancematch=0; nA=split($4,A,","); if(NR ==1) {print $2"\t"$3"\t"$4"\t"$4"\t"$4 } else{ if(nA==1){ if($1==$4){X=1} else {X=0} } else {X=0; for(a=1;a<=nA;++a) { if(A[a]==$1){X=dominance; isdominancematch=1; break}  } }; print $1"\t"$2"\t"$3"\t"X"\t"isdominancematch }  }' | grep -vw "NA$" | sed "s:$:\t$id:g" | sed '1d' > "$tmpDir"/back_"$back"_"$x".txt 
         cat "$tmpDir"/back_"$back"_"$x".txt | awk '{X+=$4}END{if(NR==0){print "00.0;0"} else {printf "%0.2f;%s\n", X/NR,NR}}' | sed "s:^:$(cat "$tmpDir"/back_"$back".txt | cut -f $x | head -1)\t:g" | sed "s:^:$x\t:g" > "$tmpDir"/back_"$back"_"$x".txt2
         ) &
         if [[ "$tcount" == $t ]]; then  wait; tcount=1; else tcount=$((tcount+1)); fi
      done
      wait
      cat "$tmpDir"/back_"$back"_*.txt2 | sort -k1,1n | cut -f 2- | sed "1iID\t$back" > "$tmpDir"/back_"$back".txt2
      if [[ "$graph" == 1 ]]
      then
         cat "$tmpDir"/back_"$back"_*.txt  | awk '{printf "%s\t%s\t%s\t%s\n", $6,$2,$3,$4}' | sed "s:$:\t$back:g" | cat <(cat "$tmpDir"/ForeGraph.graph) - > "$tmpDir"/Back."$back".graph
      fi
      )
   done
   cat "$tmpDir"/asd.txt7 | cut -f 1 > "$tmpDir"/IND.ids
   tcount=1
   for back in $(echo $BACK | tr ',' '\n' )
   do
      (
      merge "$tmpDir"/IND.ids "$tmpDir"/back_"$back".txt2 | cut -f 2- > "$tmpDir"/back_"$back".txt4
      if [[ "$graph" == 1 ]]
      then
         merge "$tmpDir"/IND.ids "$tmpDir"/Back."$back".graph | awk '{if(NR==1){print $0; X=0} else{if($1!=id){X+=1} ; printf "%04d_%s\n", X, $0 }; id=$1  }' > "$tmpDir"/Back."$back".graph2
         Rscript -e 'args = commandArgs(trailingOnly=TRUE); library(ggplot2)
         input=args[1]; output=args[2]
         A= as.data.frame(data.table::fread(input,header=TRUE))
         A$allele=as.factor(A$allele); A$qtl=as.factor(A$qtl)
         pdf(paste0(output,".pdf"),width=48,height=12)
         for (id in unique(A$ID))
         {
            a = A[which(as.character(A$ID) == id),]
            g1 <- ggplot(a,aes(pos,chr)) + geom_point(aes(size=allele,color=qtl)) + ggtitle(id) + theme(text = element_text(size=rel(6))) + theme(legend.text = element_text(size=30))
            print (g1)
         }
         dev.off()' "$tmpDir"/Back."$back".graph2 "$OutFile"_Back_"$back".pdf 2> /dev/null
      fi
      ) &
      if [[ "$tcount" == $t ]]; then  wait; tcount=1; else tcount=$((tcount+1)); fi
   done
   wait
   paste "$tmpDir"/back_*.txt4 | tr '\t' ',' | paste "$tmpDir"/asd.txt7 - > $OutFile
fi
####################################################################################
sed -i "1i$(cat $bed | cut -f 4 | cut -d',' -f 2 | tr '\n' ',' | sed "s:,$:\n:g" | sed "s:^:## :g" )" $OutFile
if [[ $clean == 1 ]]; then rm -r $tmpDir ; fi 
####################################################################################
####################################################################################
