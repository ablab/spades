/***************************************************************************
 * Title:          StatUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _STAT_UTILS_H_
#define _STAT_UTILS_H_

template<typename T>
hist(T* data, ssize_t count, ssize_t *&bins, ssize_t nbins) {
  bins = new ine[nbins];
  ssize_t i;
  for (i = 0; i < nbins; i++) {
    bins[i] = 0;
  }
  T min;
  T max;
  min = 99999999;
  max = -min;
  // Find shift of data
  for (i = 0; i < count; i++) {
    if (min > data[i])
      min = data[i];
    if (max < data[i])
      max = data[i];
  }
  
  // Bin data
  ssize_t index;
  for (i = 0; i < count; i++) {
    index = std::min(std::floor(((data[i] - min) / double(max))* nbins), nbins-1);
    bins[index]++;
  }
}
    
    
	       
  
  
  

}

#endif
