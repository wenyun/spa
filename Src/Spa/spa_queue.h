#pragma once

#include <cstdlib>
#include <functional>

class WorkQueue {
public:
  
  /// Ctor
  WorkQueue(size_t numWorkers = 0);
  
  /// Dtor, make sure there are no
  /// un-executed jobs
  ~WorkQueue();
  
  /// Dispatch a job to this queue
  void                            dispatch(std::function<void()> job);
  
  /// Partition an array and execute job for
  /// each partition
  void                            dispatch(size_t arraySize,
                                           std::function<void(size_t,size_t)> job);
  
  /// Getters
  size_t                          numWorkers() const;
  
  /// Lookup the number of cores on this compy
  static size_t                   numCores();
private:
  /// Hide implementation to avoid boost headers
  class Imp;
  
  std::unique_ptr<Imp>            _implementation;
};