#include "Lascoux.hpp"
#include "Bruhat.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

namespace {

/**
 * @brief Node in Lascoux's tree representation of the vexillary permutation partition.
 */
struct TreeNode {
    int id = 0;
    std::vector<int> children;
    int lambda_val = 0;
    int capacity = 0;

    TreeNode(int id_ = 0) : id(id_) {}
};

/**
 * @brief Constructs the Lascoux tree from permutations (x, w).
 */
std::vector<TreeNode> build_tree(const Permutation& x_perm, const Permutation& w_perm) {
    int n = x_perm.size();
    if (n != w_perm.size()) {
        throw std::invalid_argument("Lascoux: permutation size mismatch");
    }

    // Compute Lehmer code of x
    std::vector<int> code(n, 0);
    for (int i = 0; i < n; ++i) {
        int count = 0;
        for (int j = i + 1; j < n; ++j) {
            if (x_perm[j] < x_perm[i]) ++count;
        }
        code[i] = count;
    }

    std::vector<int> lambda;
    for (int v : code) {
        if (v > 0) lambda.push_back(v);
    }
    
    if (lambda.empty()) {
        return {TreeNode(0)};
    }
    std::sort(lambda.begin(), lambda.end());

    // Construct matched parenthesis string from partition lambda
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

    // Match parentheses to build the tree hierarchy
    int total_len = static_cast<int>(paren.size());
    std::vector<int> open_stack;
    std::vector<char> matched_open(total_len, 0);
    std::vector<char> matched_close(total_len, 0);

    for (int pos = 0; pos < total_len; ++pos) {
        if (paren[pos] == '(') {
            open_stack.push_back(pos);
        } else if (paren[pos] == ')' && !open_stack.empty()) {
            matched_open[open_stack.back()] = 1;
            open_stack.pop_back();
            matched_close[pos] = 1;
        }
    }

    std::vector<TreeNode> nodes;
    nodes.emplace_back(0); // Root node ID = 0
    int next_id = 1;

    std::vector<int> node_stack;
    for (int pos = 0; pos < total_len; ++pos) {
        if (paren[pos] == '(' && matched_open[pos]) {
            nodes.emplace_back(next_id);
            nodes[next_id].lambda_val = corr_lambda[pos];
            int parent = node_stack.empty() ? 0 : node_stack.back();
            nodes[parent].children.push_back(next_id);
            node_stack.push_back(next_id);
            next_id++;
        } else if (paren[pos] == ')' && matched_close[pos] && !node_stack.empty()) {
            node_stack.pop_back();
        }
    }

    // Find vertical and horizontal segments of the permutation diagram
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

    std::vector<std::pair<int, int>> crossings(n + 1);
    for (const auto& seg : h_segments) {
        int i = seg.first;
        int k = seg.second;
        for (int j : v_segments) {
            if (i <= j && j < k && x_perm[j + 1] <= x_perm[k] && x_perm[i] <= x_perm[j]) {
                int min_code = std::min(code[j], code[i]);
                if (min_code > 0) crossings[min_code] = {j, x_perm[k]};
            }
        }
    }

    // Compute rank matrix difference rxw[i][j] = rank(x[1..i, 1..j]) - rank(w[1..i, 1..j])
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

    // Assign leaf capacities from crossing rank values
    for (auto& nd : nodes) {
        if (nd.id != 0 && nd.children.empty()) {
            int c = nd.lambda_val;
            if (c > 0 && static_cast<size_t>(c) < crossings.size()) {
                auto tup = crossings[c];
                nd.capacity = (tup.first != 0 || tup.second != 0) ? rxw[tup.first][tup.second] : 0;
            } else {
                nd.capacity = 0;
            }
        }
    }

    return nodes;
}

/**
 * @brief Dynamic programming over tree to compute the generating function of edge labelings.
 */
Polynomial count_edge_labeling(const std::vector<TreeNode>& nodes) {
    int num_nodes = static_cast<int>(nodes.size());
    if (num_nodes <= 1) return Polynomial({1});

    int max_leaf_capacity = 0;
    bool found_leaf = false;
    for (const auto& node : nodes) {
        if (node.id != 0 && node.children.empty()) {
            found_leaf = true;
            max_leaf_capacity = std::max(max_leaf_capacity, node.capacity);
        }
    }
    if (!found_leaf) return Polynomial({1});

    int max_cap = max_leaf_capacity;
    std::vector<std::vector<Polynomial>> dp(num_nodes, std::vector<Polynomial>(max_cap + 2));

    auto dfs = [&](auto self, int u) -> void {
        for (int v : nodes[u].children) self(self, v);

        if (nodes[u].children.empty()) {
            int cap = std::min(max_cap + 1, nodes[u].capacity);
            for (int min_val = 0; min_val <= cap; ++min_val) {
                dp[u][min_val] = Polynomial({1});
            }
            return;
        }

        int num_children = static_cast<int>(nodes[u].children.size());
        std::vector<std::vector<Polynomial>> suffix_sums(num_children, std::vector<Polynomial>(max_cap + 2));

        for (int c_idx = 0; c_idx < num_children; ++c_idx) {
            int child = nodes[u].children[c_idx];
            Polynomial running_sum;
            for (int t = max_cap; t >= 0; --t) {
                const Polynomial& child_dp = dp[child][t];
                if (!child_dp.is_zero()) {
                    running_sum.add_shifted(child_dp, t);
                }
                suffix_sums[c_idx][t] = running_sum;
            }
        }

        for (int min_val = 0; min_val <= max_cap + 1; ++min_val) {
            Polynomial ways({1});
            bool impossible = false;

            for (int c_idx = 0; c_idx < num_children; ++c_idx) {
                const Polynomial& child_sum = suffix_sums[c_idx][min_val];
                if (child_sum.is_zero()) {
                    impossible = true;
                    break;
                }
                ways *= child_sum;
            }

            if (!impossible) {
                dp[u][min_val] = ways;
            }
        }
    };

    dfs(dfs, 0);
    return dp[0][0];
}

} // anonymous namespace

namespace Lascoux {

bool has_3412_pattern(const Permutation& p) {
    int n = p.size();
    if (n < 4) return false;

    int left_max_lt[Permutation::MAX_SIZE];
    for (int j = 0; j < n; ++j) {
        int mx = -1;
        int wj = p[j];
        for (int i = 0; i < j; ++i) {
            int wi = p[i];
            if (wi < wj && wi > mx) mx = wi;
        }
        left_max_lt[j] = mx;
    }

    int freq[Permutation::MAX_SIZE + 1] = {0};
    int pref[Permutation::MAX_SIZE + 1] = {0};

    for (int k = n - 2; k >= 1; --k) {
        ++freq[p[k + 1]];

        pref[0] = 0;
        for (int v = 1; v <= n; ++v) pref[v] = pref[v - 1] + freq[v];

        int wk = p[k];
        for (int j = 1; j < k; ++j) {
            int max_val = left_max_lt[j];
            if (max_val < 0 || max_val <= wk + 1) continue;
            if (pref[max_val - 1] - pref[wk] > 0) return true;
        }
    }

    return false;
}

Polynomial calculate_kl(const Permutation& x, const Permutation& w) {
    if (x.size() != w.size())
        throw std::invalid_argument("Lascoux: permutation size mismatch");
    if (!Bruhat::is_bruhat_less_than(x, w)) return Polynomial();
    if (x == w) {
        return Polynomial({1});
    }
    
    if (has_3412_pattern(w)) {
        throw std::invalid_argument("Lascoux method requires a 3412-avoiding permutation");
    }

    int n = x.size();
    Permutation::ValueType rev_x[Permutation::MAX_SIZE];
    Permutation::ValueType rev_w[Permutation::MAX_SIZE];
    for (int i = 0; i < n; ++i) {
        rev_x[i] = static_cast<Permutation::ValueType>(x[n - 1 - i]);
        rev_w[i] = static_cast<Permutation::ValueType>(w[n - 1 - i]);
    }

    Permutation rev_x_perm(rev_x, n);
    Permutation rev_w_perm(rev_w, n);

    auto tree = build_tree(rev_w_perm, rev_x_perm);
    return count_edge_labeling(tree);
}

} // namespace Lascoux
