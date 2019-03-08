#include <iostream>  // for cout
#include <string>  // for string
#include <algorithm>  // for replace()
#define N 500

using namespace std;

int main(void) {
    string s = "AB", temp;

    while (s.length() < N) {
        cout << s << endl;

        temp = s;
	
	    // 注意不能直接交换A, B；
        // 否则当用A代替B后，字符串中的所有字符都会变成A，这样就不能把原来的A替换成B（类似变量交换），C相当于是一个临时变量
        replace(temp.begin(), temp.end(), 'A', 'C');
        replace(temp.begin(), temp.end(), 'B', 'A');
        replace(temp.begin(), temp.end(), 'C', 'B');

        s += temp;
    }
	system("PAUSE");
    return 0;
}
