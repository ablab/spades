#include "EMPeakDetector.hh"

int EMPeakDetector::detectCoveragePeaks( MetaGraph* metaGraph, double expCovMulti[EM_PEAK_DETECTOR_MAX_NUM_CLASS], 
					 double boundarylist[EM_PEAK_DETECTOR_MAX_NUM_CLASS] ){
  double histo[EM_PEAK_DETECTOR_MAX_HISTO_LENGTH];
  double newhisto[EM_PEAK_DETECTOR_MAX_HISTO_LENGTH];
  double smoothhisto[EM_PEAK_DETECTOR_MAX_HISTO_LENGTH];
  Node *node;
  int index,secondindex = 0;
  int bin = 0;
  int loop = 0;
  int xMin = 0;
  int xMax = 1000;
  double binWidth = 1.0;
  double widthMovAve = 5.0;
  double listPeakPands[2] = {-1,-1};
  int N50 = 0;
  double thresHeight = 0;
  int thresConLen = 0;

  Graph* graph = metaGraph->getGraph();
  
  int nodeNum = nodeCount(graph);
  int *lengthlist = mallocOrExit(nodeNum, int);
  double *covlist = mallocOrExit(nodeNum, double);
  
  // Initialize expCovMulti[] and histo[]
  for (index = 0; index < 100; index++)
    expCovMulti[index] = -1;
  for (index = 0; index < EM_PEAK_DETECTOR_MAX_HISTO_LENGTH; index++)
    histo[index] = 0.0;
  
  // Make first histogram
  for (index = 1; index <= nodeNum; index++) {
    node = getNodeInGraph(graph, index);
    if (node == NULL || getNodeLength(node) <= 0)
      continue;
    bin = (int) floor( VUtils::getNodeDensity(node) );
    covlist[index - 1] = VUtils::getNodeDensity(node);
    if (bin >= EM_PEAK_DETECTOR_MAX_HISTO_LENGTH - 1)
      bin = EM_PEAK_DETECTOR_MAX_HISTO_LENGTH - 1;
    histo[bin] += getNodeLength(node);
    lengthlist[index - 1 ] = getNodeLength(node);
  }
  
  //Get N50
  N50 = n50(graph);
  printf("N50 = %d\n",N50);
  
  //Get first xMax
  thresConLen = N50 * 5;
  int maxCov = 0;
  int subMaxCov = 0;
  for (index = 0; index < nodeNum; index++){
    if (lengthlist[index] >= thresConLen){
      if (covlist[index]> maxCov){
	subMaxCov = maxCov;
	maxCov = covlist[index];
      }
    }
  }
  xMax = setXMax(subMaxCov,binWidth) + binWidth * 5;
  int firstxMax = xMax;
  free(lengthlist);
  free(covlist);
  printf("first xMAx = %d\n",firstxMax);
  
  int listWidth[2] = {0,0};
  
  while(1){
    //Set width and xMax
    setWidthByxMax(xMax,listWidth);
    binWidth = listWidth[0];
    widthMovAve = listWidth[1];
    xMax = setXMax(xMax,binWidth);
    
    //Make weighted and smoothed histogram
    xMin=0;
    weightedHisto(histo,newhisto,xMax,xMin,binWidth);
    smoothingHisto(newhisto,smoothhisto,xMax,xMin,binWidth,widthMovAve);
    
    xMin += binWidth * ((widthMovAve - 1) / 2);
    xMax -= binWidth * ((widthMovAve - 1) / 2);
    
    //Get thresHeight
    if (loop == 0){ 
      int dif = xMax - binWidth;
      thresHeight = smoothhisto[dif];
      printf("thresHeight = %f\n",thresHeight);
    }
    //Peak detection
    detectPeakPands(smoothhisto,xMax,xMin,binWidth,thresHeight,listPeakPands);
    
    
    //Record peak
    if (loop == 0)
      expCovMulti[0] = listPeakPands[0];
    expCovMulti[loop + 1] = listPeakPands[1];
    
    //When could'nt detect secondary peak, break
    if (listPeakPands[1] == -1)
      break;
    
    //Prepare for next peak
    listPeakPands[0] = listPeakPands[1];
    listPeakPands[1] = -1;
    xMax = listPeakPands[0];
    loop++;
  }
  
  if (loop == 0 && expCovMulti[0] == -1){
    puts("Error!!Couldn`t detect any peak");
    exit(1);
  }
  
  //Get parameter first value 
  int classNum = loop + 1;
  //int groupNum = floor(expCovMulti[0]) + 30;
  int groupNum = firstxMax + 30;
  long double prob[groupNum + 1][classNum];
  int count=0;
  long double ave[classNum],weight[classNum];
  
  /*for(secondindex = 1; secondindex <= groupNum; secondindex++){
    printf("%d,%f\n",secondindex,histo[secondindex]);
    } */
  
  for (index = 0; index < classNum; index++){
    ave[index] = (long double) expCovMulti[classNum - 1 -index];
    weight[index] = 1.0/classNum;
  }
  
  for(index = 0; index < classNum; index++){
    printf("Detected up-down Peak Coverage: %Lf\n",ave[index]);
  }
  //Repeat EM and judge parameter
  int preclassNum = classNum + 1;
  double d_parameter = 0.01;
  double n_parameter = 4.0;
  long double tmpave[classNum],tmpweight[classNum];
  
  while(preclassNum != classNum){
    long double preprob[groupNum + 1][EM_PEAK_DETECTOR_MAX_NUM_CLASS];
    count = 0;
    for(index = 0; index < classNum; index++){
      for(secondindex = 1; secondindex <= groupNum; secondindex++){
	preprob[secondindex][index] = 0.0;
      }
    }  
    for (index = 0; index < classNum; index++){
      weight[index] = 1.0/classNum;
    }
    
    count = myEM(histo,ave,weight,count,groupNum,classNum,preprob);
    
    /*for(index = 1; index <= groupNum; index++){
      printf("%d,",index);
      for(secondindex = 0; secondindex < classNum; secondindex++){
      printf("prob%d=%Lf,",secondindex,preprob[index][secondindex]);
      }
      printf("\n");
      }*/
    
    
    preclassNum = classNum;
    int tmpindex = 0;
    for(index = 0; index < classNum; index++){
      printf("Detected afterEM Peak Coverage: %Lf, weight: %Lf\n",ave[index],weight[index]);
    }
    for (index = 0; index < classNum; index++){
      if (weight[index] >= d_parameter){
	tmpave[tmpindex] = ave[index];
	tmpweight[tmpindex] = weight[index];
	tmpindex++;
      }
    }
    classNum = 0;
    index = 0;
    while (index < tmpindex - 1  ){ 
      judge_neighbor(index, index + 1, &index, &classNum, tmpave, tmpweight, ave, weight, n_parameter, tmpindex);
    }
    
    
    for(index = 0; index < classNum; index++){
      printf("Detected middle Peak Coverage: %Lf, weight: %Lf\n",ave[index],weight[index]);
    }
    printf("\n");
    
    if(preclassNum == classNum){
      for(index = 1; index <= groupNum; index++){
	for(secondindex = 0; secondindex < classNum; secondindex++){
	  prob[index][secondindex] = preprob[index][secondindex];
	}
      }
    }
  }
	  
  for(index = 0; index < classNum; index++){
    printf("Detected Peak Coverage: %Lf, weight: %Lf\n",ave[index],weight[index]);
  }
  
  
  //Detect the highest posterior probability
  long double max[groupNum + 1];
  int _class[groupNum + 1];
  
  for(index = 0; index <=groupNum; index++){
    max[index] = 0;
    _class[index] = 0;
  }
  
  for(index = 0; index < classNum; index++){
    for(secondindex = 1; secondindex <= groupNum; secondindex++){
      if(prob[secondindex][index] > max[secondindex]){
	max[secondindex] = prob[secondindex][index];
	_class[secondindex] = index;
      }
    }
  }	
  
  int error = 0;
  int boundarycount = 0;
  
  for (index = 0; index < 100; index++){
    expCovMulti[index] = -1;
    boundarylist[index] = -1;
  }
	
  for(index = groupNum; index > 1; index--){
    if (ave[_class[index]] < ave[_class[index - 1]]){
      error = 1;
    }
    else if (ave[_class[index]] > ave[_class[index - 1]]){
      boundarylist[boundarycount] = (double) index;
      boundarycount++;
    }
  }
  if (error == 1){
    puts("Faile EM algorithm");
    exit(1);
  }
  
  //boundary fine adjustment
  for(index = 0; index < boundarycount; index++){
    boundarylist[index] = fineboundary(histo,boundarylist[index]);
  }
  
  for(index = 0; index < classNum; index++){
    expCovMulti[classNum - 1 - index] = ave[index];
  }
  
  return classNum;
}


int EMPeakDetector::fineboundary(double *histo, int nowboundary)
{
  if(histo[nowboundary - 1] > histo[nowboundary] && histo[nowboundary] > histo[nowboundary + 1]){
    nowboundary = nowboundary + 1;
    fineboundary(histo,nowboundary);
  }
  else if(histo[nowboundary - 1] < histo[nowboundary] && histo[nowboundary] < histo[nowboundary + 1]){
    nowboundary = nowboundary - 1;
    fineboundary(histo,nowboundary);
  }
  else if(histo[nowboundary - 1] < histo[nowboundary] && histo[nowboundary] > histo[nowboundary + 1]){
    if(histo[nowboundary - 1] <= histo[nowboundary + 1]){
      nowboundary = nowboundary - 1;
    }
    else{
      nowboundary = nowboundary + 1;
    }
    fineboundary(histo,nowboundary);
  }
  else if (histo[nowboundary - 1] >= histo[nowboundary] && histo[nowboundary] <= histo[nowboundary + 1]){
    //printf("middle-nowboundary=%d\n",nowboundary);
    return (nowboundary);
  }
  
  //printf("last-nowboundary=%d\n",nowboundary);
  return (nowboundary);  
}

int EMPeakDetector::myEM(double *group,long double *newave,long double *newweight,int count,int groupNum,int classNum, long double prob[][EM_PEAK_DETECTOR_MAX_NUM_CLASS])
{
  int i,j,k;
  long double t1,t2,t3,t4,t5,t6,t7,tmp,ave[classNum],weight[classNum];
  
  if (count>32000){
    return(count);
  }
  else{
    (count)++;
  }

  for (i=0;i<classNum;i++)
    {
      ave[i] = newave[i];
      weight[i] = newweight[i];
    }
  
  for (t3=0,i=0;i<classNum;i++)
    {
      newweight[i] = 0;
      for (t4=t5=t6=t7=0,k=1;k<=groupNum;k++)
	{
	  for (t1=0,j=0;j<classNum;j++)
	    { 
	      //printf("t1=%Lf,t2=%Lf\n",t1,t2);	      
	      t2 = ave[j]/ave[i];
	      tmp=k*log(t2)+(ave[i]-ave[j])+log(weight[j]);
	      t1 += exp(tmp);
	    }
	  t6 = weight[i]/t1;
	  prob[k][i]=t6;
	  t5 = group[k]*t6;
	  newweight[i] += t5;
	  t4 += t5*k;
	  //printf("k=%d,wei=%Lf,goup[k]=%f,t1=%Lf,t6=%Lf,t5=%Lf,t4=%Lf,t7=%Lf,new=%Lf\n",k,weight[i],group[k],t1,t6,t5,t4,t7,newweight[i]);
	}
      
      newave[i]= t4/newweight[i];
      //printf("ave=%Lf,weight=%Lf\n",newave[i],newweight[i]);
      t3  += newweight[i];
    }
  //printf("\n");
  for (i=0;i<classNum;i++)
    newweight[i] /= t3 ;
  
  for (t1=0,i=0;i<classNum;i++)
    t1 += (ave[i]-newave[i])*(ave[i]-newave[i])+(weight[i]-newweight[i])*(weight[i]-newweight[i]);
  
  if ( sqrt(t1)>1e-3 )
    count = myEM(group,newave,newweight,count,groupNum,classNum,prob);
  
  return(count);
}


int EMPeakDetector::judge_neighbor(int small, int large, int *index, int *classNum, long double *tmpave,long double *tmpweight, long double *ave, long double *weight, double n_parameter, int tmpindex)
{
  if ((tmpave[large] - tmpave[small]) < n_parameter){
    if (tmpweight[large] >= tmpweight[small]){
      *index = *index + 1;
      if ( *index == (tmpindex - 1) ){
	ave[*classNum] = tmpave[large];
	weight[*classNum] = tmpweight[large];
	*classNum = *classNum + 1;
      }
    }
    else {
      *index = *index + 1;
      if (*index == (tmpindex - 1) ){
	ave[*classNum] = tmpave[small];
	weight[*classNum] = tmpweight[small];
	*classNum = *classNum + 1;
      }
      else{
	judge_neighbor(small, large + 1, index, classNum, tmpave, tmpweight, ave, weight, n_parameter, tmpindex);
      }
    }    
  }
  else {
    ave[*classNum] = tmpave[small];
    weight[*classNum] = tmpweight[small];
    *classNum = *classNum + 1;
    *index = *index + 1;
    if ( *index == (tmpindex - 1) ){
      ave[*classNum] = tmpave[large];
      weight[*classNum] = tmpweight[large];
      *classNum = *classNum + 1;
    }
  }
  return 0;
}


void EMPeakDetector::detectPeakPands(double *smoothhisto,int xMax,int xMin,int binWidth,double thresHeight,double *listPeakPands)
{
  int index,sindex;
  int countIncrease = 0; int thresIncrease = 3;
  int countDecrease = 0; int thresDecrease = 3;
  int beforeHeight = -1;
  int flagPeakStart = 0;
  int peakHeight = 0; int peakCov = 0;

  for (index = xMax - binWidth; index > xMin - binWidth; index = index - binWidth){
    if(beforeHeight == -1){
      beforeHeight = smoothhisto[index];
      continue;
    }

    if (flagPeakStart == 0){
      if (smoothhisto[index] >= thresHeight){
	if (smoothhisto[index] >= beforeHeight){
	  countIncrease++;
	  if (countIncrease >= thresIncrease){
	    countIncrease = 0;
	    flagPeakStart = 1;
	  }
	}
      }
      beforeHeight = smoothhisto[index];
    }

    else {
      if (smoothhisto[index] >= peakHeight){
	peakHeight = smoothhisto[index];
	peakCov = index;
      }
      else {
	countDecrease++;
	if (countDecrease >= thresDecrease){
	  for (sindex = 0; sindex < 2; sindex++){
	    if (listPeakPands[sindex] == -1){ 
	      listPeakPands[sindex] = peakCov + ((double) binWidth / 2 );
	      peakHeight = 0; peakCov = 0;
	      break;
	    }
	  }
	  if (listPeakPands[1] != -1) 
	    return;
	  countDecrease = 0;
	  flagPeakStart = 0;
	}
      }
    }
  }
  return;
}


void EMPeakDetector::smoothingHisto(double *newhisto,double *smoothhisto,int xMax,int xMin,int binWidth,int widthMovAve)
{
  int index,bindex;
  double binsum = 0;

  for (index = 0; index < EM_PEAK_DETECTOR_MAX_HISTO_LENGTH ; index++){
    smoothhisto[index] = 0;
  }
  
  for (index = xMin + binWidth * (widthMovAve - 1); index < xMax + binWidth; index = index + binWidth){
    for(bindex = index; bindex >= index - binWidth * (widthMovAve - 1); bindex = bindex - binWidth){
      binsum += newhisto[bindex];
    } 
    smoothhisto[index - binWidth * ((widthMovAve - 1) / 2)] = binsum / (double) widthMovAve;
    binsum = 0;
  }
}


void EMPeakDetector::weightedHisto(double *histo,double *newhisto,int xMax,int xMin,int binWidth)
{
  int index,bindex;
  for (index = 0; index < EM_PEAK_DETECTOR_MAX_HISTO_LENGTH ; index++){
    newhisto[index] = 0;
  }
  
  for (index = xMin; index < xMax + binWidth; index = index + binWidth){
    for (bindex = 0; bindex < binWidth; bindex++){
      newhisto[index] += histo[index + bindex];
    }
  }    
}


void EMPeakDetector::setWidthByxMax(int xMax,int *listWidth)
{
  if (xMax > 300){
    listWidth[0] = 6;
    listWidth[1] = 5;
  }
  if (xMax <= 300){
    listWidth[0] = 4;
    listWidth[1] = 3;
  }
  if (xMax <= 120){
    listWidth[0] = 2;
    listWidth[1] = 3;
  }
  if (xMax <= 100){
    listWidth[0] = 1;
    listWidth[1] = 1;
  }
}


int EMPeakDetector::setXMax(int xMax, int binWidth)
{
  return ((floor(xMax / binWidth)) * binWidth);
}   
