#include "spa_queue.h"

#include <mutex>

#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

using namespace std;

class WorkQueue::Imp {
public:
  
  Imp(size_t numWorkers) : _work(_ioService) {
    
    _numWorkers = numWorkers > 0 ? numWorkers : WorkQueue::numCores();
    
    for (std::size_t i = 0; i < _numWorkers; ++i) {
      _threads.create_thread(boost::bind(&boost::asio::io_service::run, &_ioService));
    }
  }
  
  ~Imp() {
    
    for (;_ioService.poll();) { }
    _ioService.stop();
    _threads.join_all();
  }
  
  void dispatch(std::function<void()> job) {
    _ioService.post(std::bind(job));
  }
  
  
  void dispatch(size_t arraySize, std::function<void(size_t,size_t)> job) {
    
    const auto numWorkers = std::min(_numWorkers, arraySize);
    
    const auto blockSize = arraySize / numWorkers;
    for (size_t worker=0;worker<numWorkers;++worker) {
      const auto startIndex = worker * blockSize;
      const auto endIndex   = (worker == (numWorkers-1)) ? arraySize : (startIndex + blockSize);
      
      _ioService.post(std::bind(job, startIndex, endIndex));
    }
  }
  
  boost::asio::io_service         _ioService;     ///< The io service managing the tasks
  boost::thread_group             _threads;       ///< The thread pool serving the io service
  boost::asio::io_service::work   _work;          ///< The work
  size_t                          _numWorkers;    ///< The number of workers we have
};

WorkQueue::WorkQueue(size_t numWorkers) {
  _implementation.reset(new Imp(numWorkers));
}

WorkQueue::~WorkQueue() {
}

void WorkQueue::dispatch(std::function<void()> job) {
  _implementation->dispatch(job);
}


void WorkQueue::dispatch(size_t arraySize, std::function<void(size_t,size_t)> job) {
  _implementation->dispatch(arraySize, job);
}

size_t WorkQueue::numWorkers() const {
  return _implementation->_numWorkers;
}

size_t WorkQueue::numCores() {
  return boost::thread::hardware_concurrency();
}



