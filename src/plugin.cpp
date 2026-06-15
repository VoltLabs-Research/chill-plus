#include <volt/plugin/plugin_entry.h>
#include <volt/plugin/chill_plus_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = ChillPlusService;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--cutoff", "O-O neighbor cutoff in angstroms", 3.5, &S::setCutoff),
};

VOLT_SERVICE_PLUGIN("volt-chill-plus", "Chill+ Ice/Hydrate Classifier", S, bindings)
