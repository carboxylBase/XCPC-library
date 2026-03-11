#include <bits/stdc++.h>
using namespace std;

int main(){
    srand(time(0));

    int n=3;

    cout<<n<<"\n";

    for(int i=1;i<=n;i++){
        cout<<rand()%20+1<<" ";
    }
    cout<<"\n";
}