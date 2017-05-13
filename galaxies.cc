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
const int MAXDOT = MAXDIM*MAXDIM;

const int UNK = -1;

int R,C;
int RD,CD;
int gridval[MAXDIM][MAXDIM][MAXDOT];
int gridgen[MAXDIM][MAXDIM][MAXDOT];
vector<bool> validgen;
int currentgen;

int D;
int dotr[MAXDOT];
int dotc[MAXDOT];

int stat_nbranches;

int getgrid(int r, int c, int d)
{
    if (!validgen[gridgen[r][c][d]]) return UNK;
    return gridval[r][c][d];
}

bool has_barrier(int r1, int c1, int r2, int c2)
{
    // TODO: This could be a little more general, in cases where we don't know exact value for either
    // but know that the sets don't intersect.
    FOR(d,D) {
        if (getgrid(r1, c1, d) == 1 && getgrid(r2, c2, d) == 0) return true;
        if (getgrid(r1, c1, d) == 0 && getgrid(r2, c2, d) == 1) return true;
    }

    return false;
}

string pretty[2*MAXDIM-1];
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
            FOR(d,D) if (getgrid(r, c, d) == 1) ch = '*';
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
        pretty[dotr[d]][dotc[d]] = 'O';
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

bool contradiction;
bool dirty;
int steps;

void setgrid_direct(int r, int c, int d, int val)
{
    gridval[r][c][d] = val;
    gridgen[r][c][d] = currentgen;
}

void setgrid(int r, int c, int d, int val)
{
    if (val == UNK) return;
    if (getgrid(r, c, d) == val) return;
    if (getgrid(r, c, d) != UNK && getgrid(r, c, d) != val) {
        contradiction = true;
        return;
    }
    setgrid_direct(r, c, d, val);
    printf("    Set (%d,%d)@%d to %d\n", r, c, d, val);
    dirty = true;
}

const int NDIR = 4;
const int DR[] = { 1, 0, -1, 0 };
const int DC[] = { 0, 1, 0, -1 };

int queue_r[MAXDIM*MAXDIM], queue_c[MAXDIM*MAXDIM];
int queue_front, queue_back;
bool mark[MAXDIM][MAXDIM];

bool enqueue(int r, int c)
{
    if (mark[r][c]) return false;
    mark[r][c] = true;
    queue_r[queue_back] = r;
    queue_c[queue_back] = c;
    ++queue_back;
    return true;
}

void prune_disconnected(int d)
{
    FOR(r,R) FOR(c,C) mark[r][c] = false;
    queue_front = queue_back = 0;

    int rd = dotr[d];
    int cd = dotc[d];

    FR(rg, rd/2, (rd+1)/2+1) {
        FR(cg, cd/2, (cd+1)/2+1) {
            enqueue(rg, cg);
        }
    }

    while (queue_front < queue_back) {
        int r = queue_r[queue_front];
        int c = queue_c[queue_front];
        ++queue_front;

        FOR(dir,NDIR) {
            int dr = DR[dir];
            int dc = DC[dir];

            int r2 = r+dr;
            int c2 = c+dc;

            if (in_bounds(r2, c2) && getgrid(r2, c2, d) != 0) {
                enqueue(r2, c2);
            }
        }
    }

    FOR(r,R) FOR(c,C) if (!mark[r][c]) setgrid(r, c, d, 0);
}

void apply_escape(int d)
{
    FOR(r,R) FOR(c,C) mark[r][c] = false;

    int nset = 0;
    FOR(r,R) FOR(c,C) if (getgrid(r, c, d) == 1) ++nset;

    bool multicomp = false;
    FOR(start_r,R) FOR(start_c,C) if (getgrid(start_r, start_c, d) == 1 && !mark[start_r][start_c]) {
        queue_front = queue_back = 0;
        enqueue(start_r, start_c);

        int escape_r = -1, escape_c = -1;
        int n_escape = 0;

        while (queue_front < queue_back) {
            int r = queue_r[queue_front];
            int c = queue_c[queue_front];
            ++queue_front;

            FOR(dir,NDIR) {
                int dr = DR[dir];
                int dc = DC[dir];

                int r2 = r+dr;
                int c2 = c+dc;

                if (in_bounds(r2, c2)) {
                    int val = getgrid(r2, c2, d);
                    if (val == 1) {
                        if (enqueue(r2, c2)) --nset;
                    } else if (val == UNK) {
                        if (n_escape == 0) {
                            escape_r = r2;
                            escape_c = c2;
                        }
                        ++n_escape;
                    }
                }
            }
        }

        if (n_escape == 1 && (nset > 0 || multicomp)) {
            setgrid(escape_r, escape_c, d, 1);
            return;
        }

        multicomp = true;
    }
}

void solve_step()
{
    printf("Solve step %d\n", steps);
    ++steps;
    dirty = false;
    contradiction = false;

    // Disjoint regions
    printf("  Disjoint regions\n");
    FOR(r,R) FOR(c,C) FOR(d,D) {
        if (getgrid(r, c, d) == 1) {
            FOR(d2,D) if (d != d2) {
                setgrid(r, c, d2, 0);
            }
        }
    }

    // Regions cover
    printf("  Regions cover\n");
    FOR(r,R) FOR(c,C) FOR(d,D) {
        int any = false;
        FOR(d2,D) if (d2 != d && getgrid(r, c, d2) != 0) any = true;
        if (!any) setgrid(r, c, d, 1);
    }

    // Rotationally symmetric regions
    printf("  Rotationally symmetric regions\n");
    FOR(r,R) FOR(c,C) FOR(d,D) {
        int r2 = r, c2 = c;
        rotate_about(dotr[d], dotc[d], &r2, &c2);

        int val = 0;
        if (in_bounds(r2, c2)) val = getgrid(r2, c2, d);
        setgrid(r, c, d, val);
    }

    // Connected regions
    printf("  Connected regions\n");
    FOR(d,D) {
        prune_disconnected(d);
    }

    // Escape
    printf("  Escape\n");
    FOR(d,D) {
        apply_escape(d);
    }

    printf("  dirty = %d\n", int(dirty));
}

bool solve()
{
    assert(validgen[currentgen]);

    steps = 0;
    do {
        solve_step();

        if (contradiction) {
            // An ancestor call made an invalid guess. Parent invalidates our generation.
            printf("Encountered contradiction\n");
            return false;
        }
    } while (dirty);

    bool complete = true;
    FOR(r,R) FOR(c,C) FOR(d,D) if (getgrid(r, c, d) == UNK) complete = false;
    if (complete) {
        printf("Completed solve!\n");
        return true;
    }

    // Try branching from any cell with minimum possibilities.
    int best_unk = MAXDOT;
    int best_r = -1;
    int best_c = -1;

    FOR(r,R) FOR(c,C) {
        int n_unk = 0;
        FOR(d,D) if (getgrid(r, c, d) == UNK) ++n_unk;
        assert(n_unk != 1);

        if (n_unk >= 2 && n_unk < best_unk) {
            best_unk = n_unk;
            best_r = r;
            best_c = c;
        }
    }

    printf("Branching from cell (%d, %d) with %d possibilities...\n", best_r, best_c, best_unk);
    ++stat_nbranches;
    FOR(d,D) if (getgrid(best_r, best_c, d) == UNK) {
        validgen.push_back(true);
        ++currentgen;
        int trialgen = currentgen;
        setgrid(best_r, best_c, d, 1);

        if (solve()) {
            return true;
        } else {
            // Clear the trial set and the non-recursive solve steps from the call.
            validgen[trialgen] = false;
        }
    }

    return false;
}

void run()
{
    stat_nbranches = 0;

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

    validgen.clear();
    validgen.push_back(false);
    currentgen = 0;

    FOR(r,R) {
        FOR(c,C) {
            FOR(d,D) {
                setgrid_direct(r, c, d, 0);
            }
        }
    }

    validgen.push_back(true);
    ++currentgen;

    FOR(d,D) {
        int rd = dotr[d];
        int cd = dotc[d];

        FR(rg, rd/2, (rd+1)/2+1) {
            FR(cg, cd/2, (cd+1)/2+1) {
                setgrid_direct(rg, cg, d, 1);
            }
        }
    }

    validgen.push_back(true);
    ++currentgen;
    bool success = solve();
    printf("Solve result :  %s\n", success ? "Success!" : "Failure");
    printf("    # branches = %d\n", stat_nbranches);

    print_grid();
}

int main()
{
    run();
    return 0;
}
