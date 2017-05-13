#include <algorithm>
#include <cassert>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <vector>
using namespace std;

#define FR(i,a,b) for(int i=(a);i<(b);++i)
#define FOR(i,n) FR(i,0,n)
#define CLR(x,a) memset(x,a,sizeof(x))
#define setmin(a,b) a = min(a,b)
#define PB push_back
#define FORALL(i,v) for(typeof((v).end())i=(v).begin();i!=(v).end();++i)
#define MP make_pair
#define A first
#define B second
#define RF(i,a,b) for(int i=(a)-1;i>=(b);--i)
#define BEND(v) (v).begin(),(v).end()
#define SZ(v) int((v).size())
#define FORI(i,v) FOR(i,SZ(v))
typedef long double ld;
typedef long long ll;

const int MAXDIM = 16;
const int MAXDOT = 256;

const int UNK = -1;

int R,C;
int RD,CD;
int grid[MAXDIM][MAXDIM][MAXDOT];

int D;
int dotr[MAXDOT];
int dotc[MAXDOT];

bool has_barrier(int r1, int c1, int r2, int c2)
{
    // TODO: This could be a little more general, in cases where we don't know exact value for either
    // but know that the sets don't intersect.
    FOR(d,D) {
        if (grid[r1][c1][d] == 1 && grid[r2][c2][d] == 0) return true;
        if (grid[r1][c1][d] == 0 && grid[r2][c2][d] == 1) return true;
    }

    return false;
}

string pretty[MAXDIM];
void print_grid()
{
    FOR(r,RD) {
        FOR(c,CD) {
            pretty[r].push_back(' ');
        }
    }

    FOR(r,R) {
        FOR(c,C) {
            char ch = '.';
            FOR(d,D) if (grid[r][c][d] == 1) ch = 'a' + d;
            pretty[2*r][2*c] = ch;
        }
    }

    FOR(r,R) {
        FOR(c,C-1) {
            if (has_barrier(r, c, r, c+1)) {
                pretty[2*r][2*c+1] = '|';
            }
        }
    }

    FOR(r,R-1) {
        FOR(c,C) {
            if (has_barrier(r, c, r+1, c)) {
                pretty[2*r+1][2*c] = '-';
            }
        }
    }

    FOR(d,D) {
        pretty[dotr[d]][dotc[d]] = 'A' + d;
    }

    FOR(r,RD) {
        printf("%s\n", pretty[r].c_str());
    }
}

bool in_bounds(int r, int c)
{
    return 0 <= r && r < R && 0 <= c && c < C;
}

void rotate_about(int rd, int cd, int *r, int *c)
{
    *r = rd - *r;
    *c = cd - *c;
}

bool dirty;

void setgrid(int r, int c, int d, int val)
{
    if (val == UNK) return;
    if (grid[r][c][d] == val) return;
    if (grid[r][c][d] != UNK) assert(grid[r][c][d] == val);
    grid[r][c][d] = val;
    printf("    Set (%d,%d)@%d to %d\n", r, c, d, val);
    dirty = true;
}

void solve_step()
{
    printf("Solve step\n");
    dirty = false;

    // Disjoint regions
    printf("  Disjoint regions\n");
    FOR(r,R) FOR(c,C) FOR(d,D) {
        if (grid[r][c][d] == 1) {
            FOR(d2,D) if (d != d2) {
                setgrid(r, c, d2, 0);
            }
        }
    }

    // Regions cover
    printf("  Regions cover\n");
    FOR(r,R) FOR(c,C) FOR(d,D) {
        int any = false;
        FOR(d2,D) if (d2 != d && grid[r][c][d2] != 0) any = true;
        if (!any) setgrid(r, c, d, 1);
    }

    // Rotationally symmetric regions
    printf("  Rotationally symmetric regions\n");
    FOR(r,R) FOR(c,C) FOR(d,D) {
        int r2 = r, c2 = c;
        rotate_about(dotr[d], dotc[d], &r2, &c2);

        int val = 0;
        if (in_bounds(r2, c2)) val = grid[r2][c2][d];
        setgrid(r, c, d, val);
    }

    // Connected regions
    // TODO

    printf("  dirty = %d\n", int(dirty));
}

void solve()
{
    do {
        solve_step();
    } while (dirty);
}

void run()
{
    scanf(" %d %d", &C, &R);
    RD = 2*R-1;
    CD = 2*C-1;

    D = 0;
    FOR(r,RD) {
        FOR(c,CD) {
            char ch;
            scanf(" %c", &ch);
            if (ch == 'o') {
                dotr[D] = r;
                dotc[D] = c;
                ++D;
            }
        }
    }

    FOR(r,R) {
        FOR(c,C) {
            FOR(d,D) {
                grid[r][c][d] = UNK;
            }
        }
    }

    FOR(d,D) {
        int rd = dotr[d];
        int cd = dotc[d];

        FR(rg, rd/2, (rd+1)/2+1) {
            FR(cg, cd/2, (cd+1)/2+1) {
                grid[rg][cg][d] = 1;
            }
        }
    }

    solve();

    print_grid();
}

int main()
{
    run();
    return 0;
}
