#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>

class Graph_M {
public:
    struct Vertex {
        std::unordered_map<std::string, int> nbrs;
    };

    std::unordered_map<std::string, Vertex> vtces;

    Graph_M() {}

    int numVertex() {
        return vtces.size();
    }

    bool containsVertex(const std::string& vname) {
        return vtces.count(vname) > 0;
    }

    void addVertex(const std::string& vname) {
        Vertex vtx;
        vtces[vname] = vtx;
    }

    void removeVertex(const std::string& vname) {
        Vertex vtx = vtces[vname];
        std::vector<std::string> keys;
        keys.reserve(vtx.nbrs.size());
        for (const auto& kvp : vtx.nbrs) {
            keys.push_back(kvp.first);
        }

        for (const auto& key : keys) {
            Vertex nbrVtx = vtces[key];
            nbrVtx.nbrs.erase(vname);
        }

        vtces.erase(vname);
    }

    int numEdges() {
        int count = 0;
        for (const auto& kvp : vtces) {
            const Vertex& vtx = kvp.second;
            count += vtx.nbrs.size();
        }
        return count / 2;
    }

    bool containsEdge(const std::string& vname1, const std::string& vname2) {
        const Vertex& vtx1 = vtces[vname1];
        const Vertex& vtx2 = vtces[vname2];
        return vtx1.nbrs.count(vname2) > 0;
    }

    void addEdge(const std::string& vname1, const std::string& vname2, int value) {
        Vertex& vtx1 = vtces[vname1];
        Vertex& vtx2 = vtces[vname2];
        if (vtx1.nbrs.count(vname2) == 0) {
            vtx1.nbrs[vname2] = value;
            vtx2.nbrs[vname1] = value;
        }
    }

    void removeEdge(const std::string& vname1, const std::string& vname2) {
        Vertex& vtx1 = vtces[vname1];
        Vertex& vtx2 = vtces[vname2];
        if (vtx1.nbrs.count(vname2) > 0) {
            vtx1.nbrs.erase(vname2);
            vtx2.nbrs.erase(vname1);
        }
    }

    void display_Map() {
        std::cout << "\t Delhi Metro Map" << std::endl;
        std::cout << "\t------------------" << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;

        for (const auto& kvp : vtces) {
            const std::string& key = kvp.first;
            const Vertex& vtx = kvp.second;
            std::string str = key + " =>\n";

            for (const auto& nbr : vtx.nbrs) {
                str += "\t" + nbr.first + "\t";
                if (nbr.first.length() < 16)
                    str += "\t";
                if (nbr.first.length() < 8)
                    str += "\t";
                str += std::to_string(nbr.second) + "\n";
            }

            std::cout << str;
        }

        std::cout << "\t------------------" << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }

    void display_Stations() {
        std::cout << "\n***********************************************************************\n";

        int i = 1;
        for (const auto& kvp : vtces) {
            std::cout << i << ". " << kvp.first << std::endl;
            i++;
        }

        std::cout << "\n***********************************************************************\n";
    }

    bool hasPath(const std::string& vname1, const std::string& vname2, std::unordered_map<std::string, bool>& processed) {
        if (containsEdge(vname1, vname2)) {
            return true;
        }

        processed[vname1] = true;

        const Vertex& vtx = vtces[vname1];
        for (const auto& nbr : vtx.nbrs) {
            if (!processed.count(nbr.first) && hasPath(nbr.first, vname2, processed)) {
                return true;
            }
        }

        return false;
    }

    struct DijkstraPair {
        std::string vname;
        std::string psf;
        int cost;

        bool operator<(const DijkstraPair& other) const {
            return cost > other.cost;
        }
    };

    int dijkstra(const std::string& src, const std::string& des, bool nan) {
        int val = 0;
        std::unordered_map<std::string, DijkstraPair> map;
        std::vector<DijkstraPair> ans;
        std::vector<std::string> keys;

        for (const auto& kvp : vtces) {
            const std::string& key = kvp.first;
            DijkstraPair np;
            np.vname = key;
            np.cost = INT_MAX;

            if (key == src) {
                np.cost = 0;
                np.psf = key;
            }

            map[key] = np;
            ans.push_back(np);
            keys.push_back(key);
        }

        std::make_heap(ans.begin(), ans.end());

        while (!ans.empty()) {
            DijkstraPair rp = ans.front();
            std::pop_heap(ans.begin(), ans.end());
ans.pop_back();
            
            if (rp.vname == des) {
                val = rp.cost;
                break;
            }

            map.erase(rp.vname);

            Vertex v = vtces[rp.vname];
            for (const auto& nbr : v.nbrs) {
                if (map.count(nbr.first) > 0) {
                    int oc = map[nbr.first].cost;
                    Vertex k = vtces[rp.vname];
                    int nc;
                    if (nan)
                        nc = rp.cost + 120 + 40 * k.nbrs[nbr.first];
                    else
                        nc = rp.cost + k.nbrs[nbr.first];

                    if (nc < oc) {
                        DijkstraPair& gp = map[nbr.first];
                        gp.psf = rp.psf + nbr.first;
                        gp.cost = nc;

                        std::make_heap(ans.begin(), ans.end());
                    }
                }
            }
        }

        return val;
    }

    struct Pair {
        std::string vname;
        std::string psf;
        int min_dis;
        int min_time;
    };

    std::string Get_Minimum_Distance(const std::string& src, const std::string& dst) {
        int min = INT_MAX;
        std::string ans = "";
        std::unordered_map<std::string, bool> processed;
        std::vector<Pair> stack;

        Pair sp;
        sp.vname = src;
        sp.psf = src + "  ";
        sp.min_dis = 0;
        sp.min_time = 0;
        stack.push_back(sp);

        while (!stack.empty()) {
            Pair rp = stack.back();
            stack.pop_back();

            if (processed.count(rp.vname) > 0) {
                continue;
            }

            processed[rp.vname] = true;

            if (rp.vname == dst) {
                int temp = rp.min_dis;
                if (temp < min) {
                    ans = rp.psf;
                    min = temp;
                }
                continue;
            }

            const Vertex& rpvtx = vtces[rp.vname];
            for (const auto& nbr : rpvtx.nbrs) {
                if (processed.count(nbr.first) == 0) {
                    Pair np;
                    np.vname = nbr.first;
                    np.psf = rp.psf + nbr.first + "  ";
                    np.min_dis = rp.min_dis + rpvtx.nbrs[nbr.first];
                    stack.push_back(np);
                }
            }
        }

        ans += std::to_string(min);
        return ans;
    }

    std::string Get_Minimum_Time(const std::string& src, const std::string& dst) {
        int min = INT_MAX;
        std::string ans = "";
        std::unordered_map<std::string, bool> processed;
        std::vector<Pair> stack;

        Pair sp;
        sp.vname = src;
        sp.psf = src + "  ";
        sp.min_dis = 0;
        sp.min_time = 0;
        stack.push_back(sp);

        while (!stack.empty()) {
            Pair rp = stack.back();
            stack.pop_back();

            if (processed.count(rp.vname) > 0) {
                continue;
            }

            processed[rp.vname] = true;

            if (rp.vname == dst) {
                int temp = rp.min_time;
                if (temp < min) {
                    ans = rp.psf;
                    min = temp;
                }
                continue;
            }

            const Vertex& rpvtx = vtces[rp.vname];
            for (const auto& nbr : rpvtx.nbrs) {
                if (processed.count(nbr.first) == 0) {
                    Pair np;
                    np.vname = nbr.first;
                    np.psf = rp.psf + nbr.first + "  ";
                    np.min_dis = rp.min_dis + rpvtx.nbrs[nbr.first];
                    np.min_time = rp.min_time + 120 + 40 * rpvtx.nbrs[nbr.first];
                    stack.push_back(np);
                }
            }
        }

        double minutes = std::ceil(static_cast<double>(min) / 60);
        ans += std::to_string(minutes);
        return ans;
    }

    std::vector<std::string> get_Interchanges(const std::string& str) {
        std::vector<std::string> arr;
        std::vector<std::string> res = splitString(str, "  ");
        arr.push_back(res[0]);
        int count = 0;
        for (int i = 1; i < res.size() - 1; i++) {
            int index = res[i].find('~');
            std::string s = res[i].substr(index + 1);

            if (s.length() == 2) {
                std::string prev = res[i - 1].substr(res[i - 1].find('~') + 1);
                std::string next = res[i + 1].substr(res[i + 1].find('~') + 1);

                if (prev == next) {
                    arr.push_back(res[i]);
                } else {
                    arr.push_back(res[i] + " ==> " + res[i + 1]);
                    i++;
                    count++;
                }
            } else {
                arr.push_back(res[i]);
            }
        }
        arr.push_back(std::to_string(count));
        arr.push_back(res[res.size() - 1Here's the C++ code converted from the given Java code:

```cpp
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <climits>
#include <cmath>

class Graph_M {
public:
    struct Vertex {
        std::unordered_map<std::string, int> nbrs;
    };

    std::unordered_map<std::string, Vertex> vtces;

    int numVetex() {
        return vtces.size();
    }

    bool containsVertex(const std::string& vname) {
        return vtces.count(vname) > 0;
    }

    void addVertex(const std::string& vname) {
        Vertex vtx;
        vtces[vname] = vtx;
    }

    void removeVertex(const std::string& vname) {
        Vertex vtx = vtces[vname];
        std::vector<std::string> keys;
        for (const auto& pair : vtx.nbrs) {
            keys.push_back(pair.first);
        }

        for (const auto& key : keys) {
            Vertex& nbrVtx = vtces[key];
            nbrVtx.nbrs.erase(vname);
        }

        vtces.erase(vname);
    }

    int numEdges() {
        int count = 0;

        for (const auto& pair : vtces) {
            const Vertex& vtx = pair.second;
            count += vtx.nbrs.size();
        }

        return count / 2;
    }

    bool containsEdge(const std::string& vname1, const std::string& vname2) {
        if (vtces.count(vname1) == 0 || vtces.count(vname2) == 0 || vtces[vname1].nbrs.count(vname2) == 0) {
            return false;
        }

        return true;
    }

    void addEdge(const std::string& vname1, const std::string& vname2, int value) {
        if (vtces.count(vname1) == 0 || vtces.count(vname2) == 0 || vtces[vname1].nbrs.count(vname2) > 0) {
            return;
        }

        vtces[vname1].nbrs[vname2] = value;
        vtces[vname2].nbrs[vname1] = value;
    }

    void removeEdge(const std::string& vname1, const std::string& vname2) {
        if (vtces.count(vname1) == 0 || vtces.count(vname2) == 0 || vtces[vname1].nbrs.count(vname2) == 0) {
            return;
        }

        vtces[vname1].nbrs.erase(vname2);
        vtces[vname2].nbrs.erase(vname1);
    }

    void display_Map() {
        std::cout << "\t Delhi Metro Map" << std::endl;
        std::cout << "\t------------------" << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;

        for (const auto& pair : vtces) {
            const std::string& key = pair.first;
            std::string str = key + " =>\n";
            const Vertex& vtx = pair.second;

            for (const auto& nbr : vtx.nbrs) {
                str += "\t" + nbr.first + "\t";
                if (nbr.first.length() < 16)
                    str += "\t";
                if (nbr.first.length() < 8)
                    str += "\t";
                str += std::to_string(nbr.second) + "\n";
            }

            std::cout << str << std::endl;
        }

        std::cout << "\t------------------" << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }

    void display_Stations() {
        std::cout << "\n***********************************************************************\n" << std::endl;
        int i = 1;

        for (const auto& pair : vtces) {
            const std::string& key = pair.first;
            std::cout << i << ". " << key << std::endl;
            i++;
        }

        std::cout << "\n***********************************************************************\n" << std::endl;
    }

    bool hasPath(const std::string& vname1, const std::string& vname2, std::unordered_map<std::string, bool>& processed) {
        if (containsEdge(vname1, vname2)) {
            return true;
        }

        processed[vname1] = true;

        Vertex vtx = vtces[vname1];
        for (const auto& nbr : vtx.nbrs) {
            if (!processed.count(nbr.first)) {
                if (hasPath(nbr.first, vname2, processed)) {
                    return true;
                }
            }
        }

        return false;
    }

    struct DijkstraPair {
        std::string vname;
        std::string psf;
        int cost;

        bool operator<(const DijkstraPair& other) const {
            return other.cost < cost;
        }
    };

    int dijkstra(const std::string& src, const std::string& des, bool nan) {
        int val = 0;
        std::vector<std::string> ans;
        std::unordered_map<std::string, DijkstraPair> map;

        std::vector<DijkstraPair> heap;

        for (const auto& pair : vtces) {
            const std::string& key = pair.first;
            DijkstraPair np;
            np.vname = key;
            np.cost = INT_MAX;

            if (key == src) {
                np.cost = 0;
                np.psf = key;
            }

            heap.push_back(np);
            map[key] = np;
        }

        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end());
            DijkstraPair rp = heap.back();
            heap.pop_back();

            if (map[rp.vname].cost != rp.cost) {
                continue;
            }

            if (rp.vname == des) {
                val = rp.cost;
                break;
            }

            Vertex vtx = vtces[rp.vname];
            for (const auto& nbr : vtx.nbrs) {
                int oc = map[nbr.first].cost;
                int nc;
                if (nan)
                    nc = rp.cost + 120 + 40 * nbr.second;
                else
                    nc = rp.cost + nbr.second;

                if (nc < oc) {
                    DijkstraPair& gp = map[nbr.first];
                    gp.vname = rp.vname;
                    gp.psf = rp.psf + nbr.first;
                    gp.cost = nc;

                    std::push_heap(heap.begin(), heap.end());
                }
            }
        }

        return val;
    }

    std::string get_Minimum_Distance(const std::string& src, const std::string& dst) {
        int min = INT_MAX;
        std::string ans = "";
        std::unordered_map<std::string, bool> processed;
        std::vector<std::string> stack;

        stack.push_back(src + "~" + std::to_string(0));
        while (!stack.empty()) {
            std::string rvtx = stack.back();
            stack.pop_back();

            if (processed.count(rvtx) > 0) {
                continue;
            }

            processed[rvtx] = true;

            std::string vname = rvtx.substr(0, rvtx.find('~'));
            int min_dis = std::stoi(rvtx.substr(rvtx.find('~') + 1));

            if (vname == dst) {
                if (min_dis < min) {
                    min = min_dis;
                    ans = rvtx;
                }
                continue;
            }

            Vertex vtx = vtces[vname];
            for (const auto& nbr : vtx.nbrs) {
                if (processed.count(nbr.first) == 0) {
                    int new_dis = min_dis + nbr.second;
                    stack.push_back(nbr.first + "~" + std::to_string(new_dis));
                }
            }
        }

        ans = ans.substr(0, ans.find('~')) + "   " + ans.substr(ans.find('~') + 1);
        ans += std::to_string(min);

        return ans;
    }

    std::string get_Minimum_Time(const std::string& src, const std::string& dst) {
        int min = INT_MAX;
        std::string ans = "";
        std::unordered_map<std::string, bool> processed;
        std::vector<std::string> stack;

        stack.push_back(src + "~" + std::to_string(0) + "~" + std::to_string(0));
        while (!stack.empty()) {
            std::string rvtx = stack.back();
            stack.pop_back();

            if (processed.count(rvtx) > 0) {
                continue;
            }

            processed[rvtx] = true;

            std::string vname = rvtx.substr(0, rvtx.find('~'));
            int min_dis = std::stoi(rvtx.substr(rvtx.find('~') + 1, rvtx.rfind('~') - rvtx.find('~') - 1));
            int min_time = std::stoi(rvtx.substr(rvtx.rfind('~') + 1));

            if (vname == dst) {
                if (min_time < min) {
                    min = min_time;
                    ans = rvtx;
                }
                continue;
            }

            Vertex vtx = vtces[vname];
            for (const auto& nbr : vtx.nbrs) {
                if (processed.count(nbr.first) == 0) {
                    int new_dis = min_dis + nbr.second;
                    int new_time = min_time + 120 + 40 * nbr.second;
                    stack.push_back(nbr.first + "~" + std::to_string(new_dis) + "~" + std::to_string(new_time));
                }
            }
        }

        double minutes = std::ceil(static_cast<double>(min) / 60);
        ans = ans.substr(0, ans.find('~')) + "   " + ans.substr(ans.find('~') + 1, ans.rfind('~') - ans.find('~') - 1);
        ans += std::to_string(minutes);

        return ans;
    }

    std::vector<std::string> get_Interchanges(const std::string& str) {
        std::vector<std::string> arr;
        std::vector<std::string> res = splitString(str, "  ");
        arr.push_back(res[0]);
        int count = 0;
        for (int i = 1; i < res.size() - 1; i++) {
            int index = res[i].find('~');
            std::string s = res[i].substr(index + 1);

            if (s.length() == 2) {
                std::string prev = res[i - 1].substr(res[i - 1].find('~') + 1);
                std::string next = res[i + 1].substr(res[i + 1].find('~') + 1);

                if (prev == next) {
                    arr.push_back(res[i]);
                } else {
                    arr.push_back(res[i] + " ==> " + res[i + 1]);
                    i++;
                    count++;
                }
            } else {
                arr.push_back(res[i]);
            }
        }
        arr.push_back(std::to_string(count));
        arr.push_back(res[res.size() - 1]);

        return arr;
    }

    std::vector<std::string> splitString(const std::string& str, const std::string& delimiter) {
        std::vector<std::string> tokens;
        std::size_t start = 0;
        std::size_t end = str.find(delimiter);

        while (end != std::string::npos) {
            tokens.push_back(str.substr(start, end - start));
            start = end + delimiter.length();
            end = str.find(delimiter, start);
        }

        tokens.push_back(str.substr(start, end));

        return tokens;
    }
};

int main() {
    Graph_M graph;

    graph.addVertex("New Delhi");
    graph.addVertex("Chandni Chowk");
    graph.addVertex("Red Fort");
    graph.addVertex("Rajiv Chowk");
    graph.addVertex("Karol Bagh");
    graph.addVertex("Paharganj");
    graph.addVertex("India Gate");
    graph.addVertex("Jama Masjid");
    graph.addVertex("Jantar Mantar");
    graph.addVertex("Lal Quila");

    graph.addEdge("New Delhi", "Chandni Chowk", 4);
    graph.addEdge("Chandni Chowk", "Red Fort", 1);
    graph.addEdge("Red Fort", "Rajiv Chowk", 3);
    graph.addEdge("Rajiv Chowk", "Karol Bagh", 3);
    graph.addEdge("Rajiv Chowk", "Paharganj", 2);
    graph.addEdgeHere's an example usage of the C++ code for the Delhi Metro Map:

```cpp
int main() {
    Graph_M graph;

    graph.addVertex("New Delhi");
    graph.addVertex("Chandni Chowk");
    graph.addVertex("Red Fort");
    graph.addVertex("Rajiv Chowk");
    graph.addVertex("Karol Bagh");
    graph.addVertex("Paharganj");
    graph.addVertex("India Gate");
    graph.addVertex("Jama Masjid");
    graph.addVertex("Jantar Mantar");
    graph.addVertex("Lal Quila");

    graph.addEdge("New Delhi", "Chandni Chowk", 4);
    graph.addEdge("Chandni Chowk", "Red Fort", 1);
    graph.addEdge("Red Fort", "Rajiv Chowk", 3);
    graph.addEdge("Rajiv Chowk", "Karol Bagh", 3);
    graph.addEdge("Rajiv Chowk", "Paharganj", 2);
    graph.addEdge("Rajiv Chowk", "India Gate", 4);
    graph.addEdge("India Gate", "Jama Masjid", 4);
    graph.addEdge("Jama Masjid", "Jantar Mantar", 2);
    graph.addEdge("Jantar Mantar", "Lal Quila", 1);
    graph.addEdge("Lal Quila", "New Delhi", 2);

    graph.display_Map();

    std::cout << "Number of stations: " << graph.numVetex() << std::endl;
    std::cout << "Number of connections: " << graph.numEdges() << std::endl;

    std::cout << "Does the metro have a path between New Delhi and Lal Quila? ";
    if (graph.hasPath("New Delhi", "Lal Quila")) {
        std::cout << "Yes" << std::endl;
    } else {
        std::cout << "No" << std::endl;
    }

    std::cout << "Minimum distance between New Delhi and Lal Quila: " << graph.get_Minimum_Distance("New Delhi", "Lal Quila") << std::endl;
    std::cout << "Minimum time between New Delhi and Lal Quila: " << graph.get_Minimum_Time("New Delhi", "Lal Quila") << std::endl;

    std::cout << "Interchanges between New Delhi and Lal Quila: ";
    std::vector<std::string> interchanges = graph.get_Interchanges(graph.get_Minimum_Distance("New Delhi", "Lal Quila"));
    for (const auto& interchange : interchanges) {
        std::cout << interchange << " ";
    }
    std::cout << std::endl;

    return 0;
}
