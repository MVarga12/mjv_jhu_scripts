/*
 * Created:  8-May-2018 02:14:57 PM EDT
 * Modified:  8-May-2018 02:40:15 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>

namespace Involvement {
/* Involved ifaces are reactant/product ifaces which directly change in the reaction
 * can be binding ifaces, unbinding ifaces, or ifaces which only change state
 */
enum Involved { binding, unbinding, onlyStateChange };

/* Ancillary ifaces are ifaces in reactant or product molecules which are not reactant or product ifaces
 * can be required bound/unbound ifaces or required states
 */
enum Ancillary { bound, unbound, state };
}

int main() {}
