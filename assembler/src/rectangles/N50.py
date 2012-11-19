############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

def NG50(numlist, reference_length, percentage = 50.0):
  """
  Abstract: Returns the NG50 value of the passed list of numbers.
  Comments: Works for any percentage (e.g. NG60, NG70) with optional argument
  Usage: NG50(numlist, reference_length)

  Based on the definition from this SEQanswers post
  http://seqanswers.com/forums/showpost.php?p=7496&postcount=4
  (modified Broad Institute's definition
  https://www.broad.harvard.edu/crd/wiki/index.php/N50)
  
  See SEQanswers threads for details:
  http://seqanswers.com/forums/showthread.php?t=2857
  http://seqanswers.com/forums/showthread.php?t=2332
  """
  assert percentage >= 0.0
  assert percentage <= 100.0
  numlist.sort(reverse = True)
  s = reference_length
  limit = reference_length * (100.0 - percentage) / 100.0
  for l in numlist:
    s -= l
    if s <= limit:
      return l


def N50(numlist, percentage = 50.0):
  """
  Abstract: Returns the N50 value of the passed list of numbers.
  Comments: Works for any percentage (e.g. N60, N70) with optional argument
  Usage: N50(numlist)
  """
  return NG50(numlist, sum(numlist), percentage)

