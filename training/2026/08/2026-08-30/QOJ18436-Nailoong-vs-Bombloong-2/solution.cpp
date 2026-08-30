#include <bits/stdc++.h>

using namespace std;

void solveFirstRun() {
    int n;
    cin >> n;

    vector<vector<int>> graph(n + 1);
    for (int i = 1; i < n; ++i) {
        int u, v;
        cin >> u >> v;
        graph[u].push_back(v);
        graph[v].push_back(u);
    }

    int root = 1;
    for (int i = 1; i <= n; ++i) {
        if (graph[i].size() == 1) {
            root = i;
            break;
        }
    }

    queue<int> que;
    vector<int> order;
    vector<bool> visited(n + 1, false);

    que.push(root);
    visited[root] = true;

    while (!que.empty()) {
        int u = que.front();
        que.pop();
        order.push_back(u);

        for (int v : graph[u]) {
            if (!visited[v]) {
                visited[v] = true;
                que.push(v);
            }
        }
    }

    cout << n - 3 << '\n';
    for (int i = 1; i <= n - 3; ++i) {
        cout << order[i] << " \n"[i == n - 3];
    }
}

void solveSecondRun() {
    int n, m;
    cin >> n >> m;

    vector<int> cut(n + 1, 0);
    for (int i = 2; i <= n - 2; ++i) {
        cin >> cut[i];
    }

    vector<vector<int>> graph(n + 1);
    int remainingDegreeSum = 2 * n - 2;

    graph[1].push_back(2);
    remainingDegreeSum -= 2;

    vector<int> currentLevel;
    for (int i = 1; i <= cut[2] - 1; ++i) {
        graph[2].push_back(i + 2);
        remainingDegreeSum -= 2;
        currentLevel.push_back(i + 2);
    }

    int lastVertex = cut[2] + 1;

    while (!currentLevel.empty()) {
        vector<int> nextLevel;

        for (int v : currentLevel) {
            if (v > n - 2) continue;

            int children = cut[v] - cut[v - 1] + 1;
            while (children--) {
                graph[v].push_back(++lastVertex);
                remainingDegreeSum -= 2;
                nextLevel.push_back(lastVertex);
            }
        }

        currentLevel = move(nextLevel);
    }

    assert(remainingDegreeSum == 0 || remainingDegreeSum == 2);
    if (remainingDegreeSum == 2) {
        graph[n - 1].push_back(n);
    }

    for (int u = 1; u <= n; ++u) {
        for (int v : graph[u]) {
            cout << u << ' ' << v << '\n';
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int run;
    cin >> run;

    if (run == 1) {
        solveFirstRun();
    } else {
        solveSecondRun();
    }

    return 0;
}
