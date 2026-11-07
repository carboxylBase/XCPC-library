#include <bits/stdc++.h>
using namespace std;

int main(){
    freopen("input.txt", "w", stdout);
    int n, q;
    n = 500000;
    q = 1e5;
    cout << 1 << endl;
    cout << n << ' ' << q << endl;
    for (int i = 2;i<n+1;i++) cout << i - 1 << ' ';
    cout << endl;
    while (q --) {
        int x = rand() % n + 1;
        cout << x << endl;
    } 
}