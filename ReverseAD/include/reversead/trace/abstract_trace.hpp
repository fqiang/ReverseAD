#ifndef ABSTRACT_TRACE_H_
#define ABSTRACT_TRACE_H_

#include "reversead/logger.hpp"
#include "reversead/reversead_base.hpp"

namespace ReverseAD {

class AbstractTrace {
 public:
  AbstractTrace() {
    num_ind = 0;
    num_dep = 0;
    num_dummy_ind = 0;
    num_dummy_dep = 0;
  }
  virtual ~AbstractTrace() {}

  inline virtual void put_op(const opbyte&) = 0;
  inline virtual void put_loc(const locint&) = 0;
  inline virtual void put_val(const double&) = 0;
  inline virtual void put_sr_info(const SendRecvInfo&) = 0;
  
  // forward sweep
  virtual void init_forward() = 0;
  virtual void end_forward() = 0;
  virtual opbyte get_next_op_f() = 0;
  virtual locint get_next_loc_f() = 0;
  virtual double get_next_val_f() = 0;
  virtual bool has_next_sr_info_f() = 0;
  virtual SendRecvInfo get_next_sr_info_f() = 0;

  // reverse sweep
  virtual void init_reverse() = 0;
  virtual void end_reverse() = 0;
  virtual opbyte get_next_op_r() = 0;
  virtual locint get_next_loc_r() = 0;
  virtual double get_next_val_r() = 0;
  virtual SendRecvInfo get_next_sr_info_r() = 0;

  // for debug
  virtual void dump_trace() = 0;
  virtual void dump_trace(Logger& logger) = 0;
  
  void declare_ind() {num_ind++;}
  void declare_dep() {num_dep++;}
  int get_num_ind() {return num_ind;}
  int get_num_dep() {return num_dep;}
  void increase_dummy_ind(int size) {num_dummy_ind += size;}
  void increase_dummy_dep(int size) {num_dummy_dep += size;}
  int get_num_dummy_ind() {return num_dummy_ind;}
  int get_num_dummy_dep() {return num_dummy_dep;}
 private:
  int num_ind;
  int num_dep;
  int num_dummy_ind;
  int num_dummy_dep;
};

} // namespace ReverseAD

#endif // ABSTRACT_TRACT_H_