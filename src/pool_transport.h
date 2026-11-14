// The worker pool's byte transport (src/pool_transport.cpp) for native
// callers: a message's payload over one of the pool's connections.
#ifndef EMC2_POOL_TRANSPORT_H
#define EMC2_POOL_TRANSPORT_H

#include <Rcpp.h>
#include <chrono>
#include <functional>
#include <string>
#include <vector>

namespace emc_pool {

// How a read waits. Untimed, an empty read means the other end closed. Timed,
// it polls up to `limit` seconds, asking `alive` now and then whether the
// other end still runs and, at the limit, `busy` whether it is still working
// (which restarts the wait). An unset function is not asked.
struct Watch {
  bool timed = false;
  double limit = 0.0;
  std::function<bool()> alive, busy;
};

// The connection behind an R connection object; NULL when it is not one.
void* connection(SEXP con);
// One message's payload from the connection, its segments joined. An empty
// string on success, else what went wrong.
std::string recv_payload(void* con, const Watch& watch, std::vector<unsigned char>& payload);
// A message's payload to the connection in segments, then flushed. An empty
// string on success. A timed write gives up after `limit` seconds.
std::string send_payload(void* con, const unsigned char* bytes, size_t n, bool timed,
                         double limit);
// Up to n bytes that are ready on a non-blocking connection, without waiting.
size_t read_ready(void* con, unsigned char* out, size_t n);
// One wait of a polling loop that has been idle since `idle_since`: nothing
// at first, so a reply about to arrive is taken at once, then a sleep that
// backs off. `wait` carries the sleep between calls (start it at zero).
void idle_wait(const std::chrono::steady_clock::time_point& idle_since, double& wait);

}  // namespace emc_pool

#endif
