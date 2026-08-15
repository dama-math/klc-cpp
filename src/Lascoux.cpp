#include "Lascoux.hpp"
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <functional>
#include <iostream>

namespace {

/**
 * @brief Internal node structure for Lascoux's tree representation.
 */
struct Node {
    std::vector<int> children; ///< List of child node IDs.
    int id;                    ///< Unique identifier for the node.
    int corr_c = 0;            ///< Corresponding lambda value from the construction phase.
    int capacity = 0;          ///< Capacity of the node (relevant for leaves).
    
    /**
     * @brief Constructor for a new Node.
     * @param id_ The ID to assign to the node.
     */
    Node(int id_ = 0) : id(id_) {}
};

/**
 * @brief Builds the tree structure used in Lascoux's algorithm from two permutations.
 * 
 * @param x_perm The first permutation (typically x).
 * @param w_perm The second permutation (typically w).
 * @return A vector of Node objects representing the built tree.
 */
std::vector<Node> build_tree_from_perms(const Permutation& x_perm, const Permutation& w_perm) {
    int n = x_perm.size();

    if (n != w_perm.size()) {
        throw std::invalid_argument("Lascoux::build_tree_from_perms: permutations of different sizes");
    }

    // Step 1: Build Lehmer code and lambda sequence from x
    std::vector<int> code(n, 0);
    for (int i = 0; i < n; ++i) {
        int cnt = 0;
        for (int j = i + 1; j < n; ++j) {
            if (x_perm[j] < x_perm[i]) ++cnt;
        }
        code[i] = cnt;
    }

    std::vector<int> lambda;
    for (int v : code) {
        if (v > 0) lambda.push_back(v);
    }
    
    // Return a single root node if there are no descents
    if (lambda.empty()) {
        return {Node(0)}; 
    }
    std::sort(lambda.begin(), lambda.end());

    // Step 2: Build the parenthesis word representation from lambda
    std::string paren;
    std::vector<int> corr_lambda;

    paren += std::string(lambda[0], '(');
    for (int i = 0; i < lambda[0]; ++i) corr_lambda.push_back(lambda[0]);
    paren += ')';
    corr_lambda.push_back(lambda[0]);
    for (size_t i = 1; i < lambda.size(); ++i) {
        int horizontal = lambda[i] - lambda[i - 1];
        paren += std::string(horizontal, '(');
        for (int j = 0; j < horizontal; ++j) corr_lambda.push_back(lambda[i]);
        paren += ')';
        corr_lambda.push_back(lambda[i]);
    }

    // Step 3: Construct the tree from the parenthesis word
    std::vector<Node> nodes;
    nodes.emplace_back(0); // Root node has ID = 0
    int next_id = 1;

    int L = (int)paren.size();
    std::vector<int> open_stack;
    std::vector<char> matched_open(L, 0);
    std::vector<char> matched_close(L, 0);

    for (int pos = 0; pos < L; ++pos) {
        if (paren[pos] == '(') {
            open_stack.push_back(pos);
        } else if (paren[pos] == ')') {
            if (!open_stack.empty()) {
                matched_open[open_stack.back()] = 1;
                open_stack.pop_back();
                matched_close[pos] = 1;
            }
        }
    }

    std::vector<int> node_stack;
    for (int pos = 0; pos < L; ++pos) {
        if (paren[pos] == '(' && matched_open[pos]) {
            nodes.emplace_back(next_id);
            nodes[next_id].corr_c = corr_lambda[pos];
            int parent = node_stack.empty() ? 0 : node_stack.back();
            nodes[parent].children.push_back(next_id);
            node_stack.push_back(next_id);
            next_id++;
        } else if (paren[pos] == ')' && matched_close[pos] && !node_stack.empty()) {
            node_stack.pop_back();
        }
    }

    // Step 4: Identify vertical and horizontal segments from x
    std::vector<int> v_segments;
    for (int i = 0; i < n - 1; ++i) {
        if (x_perm[i] > x_perm[i + 1]) v_segments.push_back(i);
    }

    std::vector<int> x_inv(n, 0);
    for (int i = 0; i < n; ++i) x_inv[x_perm[i]] = i;

    std::vector<std::pair<int, int>> h_segments;
    for (int j = 0; j < n - 1; ++j) {
        if (x_inv[j + 1] < x_inv[j]) {
            h_segments.emplace_back(x_inv[j + 1], x_inv[j]);
        }
    }

    // Step 5: Identify segment crossings
    std::vector<std::pair<int, int>> crossings(n + 1); // Uses 1-based indexing for code values
    for (const auto& seg : h_segments) {
        int i = seg.first;
        int k = seg.second;
        for (int j : v_segments) {
            if (i <= j && j < k && x_perm[j + 1] <= x_perm[k] && x_perm[i] <= x_perm[j]) {
                int min_code = std::min(code[j], code[i]);
                if(min_code > 0) crossings[min_code] = {j, x_perm[k]};
            }
        }
    }

    // Step 6: Define capacities for the leaf nodes based on the rank matrix differences
    std::vector<std::vector<int>> rxw(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        rxw[i][x_perm[i]]++;
        rxw[i][w_perm[i]]--;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 1; j < n; ++j) rxw[i][j] += rxw[i][j - 1];
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 1; i < n; ++i) rxw[i][j] += rxw[i - 1][j];
    }

    for (auto& nd : nodes) {
        if (nd.id != 0 && nd.children.empty()) { // Process leaf nodes (excluding the root)
            int c = nd.corr_c;
            if (c > 0 && static_cast<size_t>(c) < crossings.size()) {
                auto tup = crossings[c];
                if (tup.first != 0 || tup.second != 0) {
                    nd.capacity = rxw[tup.first][tup.second];
                } else {
                    nd.capacity = 0;
                }
            } else {
                 nd.capacity = 0;
            }
        }
    }

    return nodes;
}

/**
 * @brief Computes the polynomial associated with the valid edge labelings of the tree.
 * 
 * @param nodes The tree represented as a vector of Node objects.
 * @return The computed Kazhdan-Lusztig polynomial.
 */
Polynomial count_edge_labeling(const std::vector<Node>& nodes) {
    int n = nodes.size();
    if (n <= 1) return Polynomial({1});

    std::vector<int> parent(n, -1);
    int root = -1;
    for (int u = 0; u < n; ++u) {
        for (int v : nodes[u].children) parent[v] = u;
    }
    for (int i = 0; i < n; ++i) {
        if (parent[i] == -1) {
            root = i;
            break;
        }
    }
    if (root == -1) root = 0;

    int max_leaf_capacity = 0;
    bool found_leaf = false;
    for (const auto& node : nodes) {
        if (node.id != 0 && node.children.empty()) {
            found_leaf = true;
            max_leaf_capacity = std::max(max_leaf_capacity, node.capacity);
        }
    }
    if (!found_leaf) return Polynomial({1});

    int C = max_leaf_capacity;
    std::vector<std::vector<Polynomial>> dp(n, std::vector<Polynomial>(C + 2));

    // DFS to compute the polynomials bottom-up
    std::function<void(int)> dfs = 
        [&](int u) {
        for (int v : nodes[u].children) dfs(v);

        if (nodes[u].children.empty()) {
            for (int minVal = 0; minVal <= C + 1; ++minVal) {
                if (minVal <= nodes[u].capacity) {
                    dp[u][minVal] = Polynomial({1});
                } else {
                    dp[u][minVal] = Polynomial(); 
                }
            }
            return;
        }

        for (int minVal = 0; minVal <= C + 1; ++minVal) {
            Polynomial waysHere({1});
            bool impossible = false;

            for (int child : nodes[u].children) {
                Polynomial childSum;
                for (int t = minVal; t <= C; ++t) {
                    const Polynomial& childDP = dp[child][t];
                    if (childDP.is_zero()) continue;
                    childSum += childDP.shift(t);
                }
                if (childSum.is_zero()) {
                    impossible = true;
                    break;
                }
                waysHere *= childSum;
            }

            if (impossible) {
                dp[u][minVal] = Polynomial();
            } else {
                dp[u][minVal] = waysHere;
            }
        }
    };

    dfs(root);

    return dp[root][0];
}

} // anonymous namespace

namespace Lascoux {

/**
 * @brief Checks if a permutation contains a 3412 pattern.
 * 
 * A permutation has a 3412 pattern if there exist indices i < j < k < l
 * such that w[k] < w[l] < w[i] < w[j].
 * 
 * @param p The permutation to check.
 * @return True if the pattern exists, false otherwise.
 */
bool has_3412_pattern(const Permutation& p) {
    int n = p.size();
    if (n < 4) return false;

    int left_max_lt[32];
    for (int j = 0; j < n; ++j) {
        int mx = -1;
        int wj = p[j];
        for (int i = 0; i < j; ++i) {
            int wi = p[i];
            if (wi < wj && wi > mx) mx = wi;
        }
        left_max_lt[j] = mx;
    }

    int freq[33] = {0};
    int pref[33] = {0};

    for (int k = n - 2; k >= 1; --k) {
        ++freq[p[k + 1]];

        pref[0] = 0;
        for (int v = 1; v <= n; ++v) pref[v] = pref[v - 1] + freq[v];

        int wk = p[k];
        for (int j = 1; j < k; ++j) {
            int T = left_max_lt[j];
            if (T < 0 || T <= wk + 1) continue;
            // Check if there is a valid 'l' element completing the pattern
            if (pref[T - 1] - pref[wk] > 0) return true;
        }
    }

    return false;
}

/**
 * @brief Calculates the Kazhdan-Lusztig polynomial using Lascoux's algorithm.
 * 
 * @param x The inner permutation.
 * @param w The outer permutation.
 * @return The calculated Kazhdan-Lusztig polynomial.
 */
Polynomial calculate_kl(const Permutation& x, const Permutation& w) {
    if (x == w) {
        return Polynomial({1});
    }
    
    // Lascoux algorithm is only valid for 3412 pattern-avoiding permutations
    if (has_3412_pattern(w)) {
        std::cerr << "Warning: 3412 pattern detected. Lascoux method invalid.\n";
        return Polynomial(); 
    }

    // Convert to vectors and reverse to apply the tree construction
    std::vector<int> x_data = x.to_vector();
    std::vector<int> w_data = w.to_vector();
    std::reverse(x_data.begin(), x_data.end());
    std::reverse(w_data.begin(), w_data.end());

    Permutation rev_x_perm(x_data);
    Permutation rev_w_perm(w_data);

    // Build tree and count valid edge labelings
    auto tree = build_tree_from_perms(rev_w_perm, rev_x_perm);
    return count_edge_labeling(tree);
}

} // namespace Lascoux
