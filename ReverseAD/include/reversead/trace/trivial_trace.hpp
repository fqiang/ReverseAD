#ifndef TRIVIAL_TRACE_H_
#define TRIVIAL_TRACE_H_

#include <iostream>
#include <memory>

#include "reversead/common/reversead_config.h"
#include "reversead/common/reversead_type.hpp"
#include "reversead/tape/trivial_tape.hpp"
#include "reversead/tape/disk_tape.hpp"
#include "reversead/trace/abstract_trace.hpp"

namespace ReverseAD {

#ifdef ENABLE_DISK_TAPE
template <typename T>
using VirtualTape = DiskTape<T>;
#else
template <typename T>
using VirtualTape = TrivialTape<T>;
#endif

template<typename Base>
class TrivialTrace : public AbstractTrace<Base> {
  using AbstractTrace<Base>::num_ind;
  using AbstractTrace<Base>::num_dep;
  using AbstractTrace<Base>::num_dummy_ind;
  using AbstractTrace<Base>::num_dummy_dep;
  using AbstractTrace<Base>::num_param;

 public:
  template <typename OldBase, typename NewBase>
  friend std::shared_ptr<TrivialTrace<NewBase>> copy_tape(
    const std::shared_ptr<TrivialTrace<OldBase>>& other,
    const std::shared_ptr<VirtualTape<NewBase>>& val_tape,
    const std::shared_ptr<VirtualTape<NewBase>>& param_tape);
  template <typename OldBase, typename NewBase>
  friend std::shared_ptr<TrivialTrace<NewBase>> copy_tape(
    const std::shared_ptr<TrivialTrace<OldBase>>& other,
    const std::shared_ptr<VirtualTape<NewBase>>& val_tape);

  // Default c-tor and d-tor
  TrivialTrace() {
    op_tape = std::make_shared<VirtualTape<opbyte>>();
    loc_tape = std::make_shared<VirtualTape<locint>>();
    val_tape = std::make_shared<VirtualTape<Base>>();
    param_tape = std::make_shared<VirtualTape<Base>>();
    coval_tape = std::make_shared<VirtualTape<double>>();
#ifdef ENABLE_REVERSEAD_MPI
    sr_info_tape = std::make_shared<TrivialTape<SendRecvInfo>>();
    comm_loc_tape = std::make_shared<TrivialTape<locint>>();
#endif
  }
  ~TrivialTrace() {
    op_tape.reset();
    loc_tape.reset();
    val_tape.reset();
    param_tape.reset();
    coval_tape.reset();
#ifdef ENABLE_REVERSEAD_MPI
    sr_info_tape.reset();
    comm_loc_tape.reset();
#endif
  }

  void init_tracing() {
    op_tape->init_taping();
    loc_tape->init_taping();
    val_tape->init_taping();
    param_tape->init_taping();
    coval_tape->init_taping();
#ifdef ENABLE_REVERSEAD_MPI
    sr_info_tape->init_taping();
    comm_loc_tape->init_taping();
#endif
  }
  void end_tracing() {
    op_tape->end_taping();
    loc_tape->end_taping();
    val_tape->end_taping();
    param_tape->end_taping();
    coval_tape->end_taping();
#ifdef ENABLE_REVERSEAD_MPI
    sr_info_tape->end_taping();
    comm_loc_tape->end_taping();
#endif
  }
  // Write
  inline void put_op(const opbyte& opcode) {
    op_tape->put(opcode);
  }
  inline void put_loc(const locint& loc) {
    loc_tape->put(loc);
  }
  inline void put_val(const Base& val) {
    val_tape->put(val);
  }
  inline void put_param(const Base& param) {
    num_param++;
    param_tape->put(param);
  }
  inline void put_coval(const double& coval) {
    coval_tape->put(coval);
  }

  inline void put_sr_info(const SendRecvInfo& sr_info) {
    sr_info_tape->put(sr_info);
  }
  inline void put_comm_loc(const locint& comm_loc) {
    std::cout << "trace->put_comm_loc : " << comm_loc << std::endl;
    comm_loc_tape->put(comm_loc);
  }
  inline void init_comm_forward() {
    sr_info_tape->init_forward();
    comm_loc_tape->init_forward();
  }
  inline void end_comm_forward() {
    sr_info_tape->end_forward();
    comm_loc_tape->end_forward();
  }
  inline bool has_next_sr_info_f() {
    return sr_info_tape->has_next_f();
  }
  inline SendRecvInfo get_next_sr_info_f() {
    return sr_info_tape->get_next_f();
  }
  inline locint get_next_comm_loc_f() {
    return comm_loc_tape->get_next_f();
  }
 
  // forward sweep
  inline void init_forward() {
    op_tape->init_forward();
    loc_tape->init_forward();
    val_tape->init_forward();
    param_tape->init_forward();
    coval_tape->init_forward();
  }
  inline void end_forward() {
    op_tape->end_forward();
    loc_tape->end_forward();
    val_tape->end_forward();
    param_tape->end_forward();
    coval_tape->end_forward();
  }
  inline opbyte get_next_op_f() {
    return op_tape->get_next_f();
  }
  inline locint get_next_loc_f() {
    return loc_tape->get_next_f();
  }
  inline Base get_next_val_f() {
    return val_tape->get_next_f();
  }
  inline Base get_next_param_f() {
    return param_tape->get_next_f();
  }
  inline double get_next_coval_f() {
    return coval_tape->get_next_f();
  }

  // reverse sweep
  inline void init_reverse() {
    op_tape->init_reverse();
    loc_tape->init_reverse();
    val_tape->init_reverse();
    param_tape->init_reverse();
    coval_tape->init_reverse();
  }
  inline void end_reverse() {
    op_tape->end_reverse();
    loc_tape->end_reverse();
    val_tape->end_reverse();
    param_tape->end_reverse();
    coval_tape->end_reverse();
  }
  inline opbyte get_next_op_r() {
    return op_tape->get_next_r();
  }
  inline locint get_next_loc_r() {
    return loc_tape->get_next_r();
  }
  inline Base get_next_val_r() {
    return val_tape->get_next_r();
  }
  inline Base get_next_param_r() {
    return param_tape->get_next_r();
  }
  inline double get_next_coval_r() {
    return coval_tape->get_next_r();
  }

  // for debug
  inline void dump_trace() {
    std::cout << "Op tape:" << std::endl;
    op_tape->dump_tape();
    std::cout << "Loc tape:" << std::endl;
    loc_tape->dump_tape();
    std::cout << "Val tape:" << std::endl;
    val_tape->dump_tape();
    std::cout << "Param tape:" << std::endl;
    param_tape->dump_tape();
    std::cout << "Const tape:" << std::endl;
    coval_tape->dump_tape();
    if (sr_info_tape) {
      std::cout << "SendRecvInfo tape:" << std::endl;
      sr_info_tape->dump_tape();
      std::cout << "Comm Info tape:" << std::endl;
      comm_loc_tape->dump_tape();
    }
  }
   
 private:
  std::shared_ptr<VirtualTape<opbyte>> op_tape;
  std::shared_ptr<VirtualTape<locint>> loc_tape;
  std::shared_ptr<VirtualTape<Base>> val_tape;
  std::shared_ptr<VirtualTape<Base>> param_tape;
  std::shared_ptr<VirtualTape<double>> coval_tape;
  std::shared_ptr<TrivialTape<SendRecvInfo>> sr_info_tape;
  std::shared_ptr<TrivialTape<locint>> comm_loc_tape;
};

template <typename OldBase, typename NewBase>
std::shared_ptr<TrivialTrace<NewBase>> copy_tape(
      const std::shared_ptr<TrivialTrace<OldBase>>& other,
      const std::shared_ptr<VirtualTape<NewBase>>& val_tape,
      const std::shared_ptr<VirtualTape<NewBase>>& param_tape) {
    std::shared_ptr<TrivialTrace<NewBase>> ret =
        std::make_shared<TrivialTrace<NewBase>>();
    ret->op_tape = other->op_tape;
    ret->loc_tape = other->loc_tape;
    ret->val_tape = val_tape;
    ret->param_tape = param_tape;
    ret->coval_tape = other->coval_tape;
    ret->num_ind = other->num_ind;
    ret->num_dep = other->num_dep;
#ifdef ENABLE_REVERSEAD_MPI
    ret->sr_info_tape = other->sr_info_tape;
    ret->comm_loc_tape = other->comm_loc_tape;
    ret->num_dummy_ind = other->num_dummy_ind;
    ret->num_dummy_dep = other->num_dummy_dep;
#endif
    return ret;
  }

template <typename OldBase, typename NewBase>
static std::shared_ptr<TrivialTrace<NewBase>> copy_tape(
      const std::shared_ptr<TrivialTrace<OldBase>>& other,
      const std::shared_ptr<VirtualTape<NewBase>>& val_tape) {
    // TODO(some warning message when other trace has non empty parameter)
    return copy_tape<OldBase, NewBase>(
        other, val_tape, std::make_shared<VirtualTape<NewBase>>());
}

} // namespace ReverseAD

#endif // TRIVIAL_TRACE_H_
