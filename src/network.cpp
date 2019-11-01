#include "network.h"
#include "random.h"
#include <map>

std::set<size_t> Network::step(const std::vector<double>& means_){
	//------------------------------ à coder -----------------------------//

	std::set<size_t> set_;
	double w(0.0);
	double exitation(0.0);
	double inhibition(0.0);
	
	for(size_t i(0); i < neurons.size(); ++i){
			if(neurons[i].firing()) {
				neurons[i].reset();
				set_.insert(i);
			}
	}
	
	
		for(size_t i = 0 ; i < neurons.size(); ++i){
				if(neurons[i].is_inhibitory()){
					w = 0.4;
				} else {
					w = 1.0;
				}
				std::vector<std::pair<size_t, double> > neighb = neighbors(i);
				for(auto voisin : neighb){
					if (set_.set::find(voisin.first) != set_.end()){
						if(neurons[voisin.first].is_inhibitory()){
							inhibition += voisin.second;
						} else {
							exitation += voisin.second;
						}
					}
				}
			neurons[i].input(w*means_[i] + 0.5*exitation + inhibition);
			exitation = 0.0;
			inhibition = 0.0;
			neurons[i].step();
			}
	return set_;
}

std::pair<size_t, double> Network::degree(const size_t& index_) const{
	//------------------------------ à coder -----------------------------//
	std::vector<std::pair<size_t, double> > neighbors_(neighbors(index_));
	size_t nb_connections (neighbors_.size());
	double sum_intensity(0.0);
	for(const auto& sm_pair : neighbors_){
			sum_intensity += sm_pair.second;
	}
	return {nb_connections, sum_intensity};
}


std::vector<std::pair<size_t, double> > Network::neighbors(const size_t& nb) const{
	//------------------------------ à coder -----------------------------//
	std::vector<std::pair<size_t, double> > neighbors_;
	std::map<std::pair<size_t, size_t>, double>::const_iterator itlow;
	itlow = links.lower_bound ({nb, 0});
	for (auto it = itlow; it != links.end() and it->first.first == nb; ++it){
			neighbors_.push_back({it->first.second, it->second});
	}
	return neighbors_;
}

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}
