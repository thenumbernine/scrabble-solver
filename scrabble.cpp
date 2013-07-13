/*
	Copyright (c) 2010 Christopher E. Moore ( christopher.e.moore@gmail.com / http://christopheremoore.net )

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in
	all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
	THE SOFTWARE.
*/
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <memory.h>
#include <assert.h>

#ifdef WIN32
#include <windows.h>	//CreateThread and all that jazz
#include <time.h>
#define strcasecmp _stricmp
#else
#include <sys/time.h>
#endif

#define NUMBER_OF_THREADS 4
#define BOARDSIZE	15
#define MAXHANDSIZE	7

#define numberof(x)	(sizeof(x)/sizeof((x)[0]))

using namespace std;

class vec2i {
public: 
	int x,y; 
	vec2i() : x(0), y(0) {} 
	vec2i(int _x, int _y) : x(_x), y(_y) {}
	vec2i &operator=(const vec2i &v) { x = v.x; y = v.y; return *this; }
};
vec2i &operator+=(vec2i &a, const vec2i &b) { a.x += b.x; a.y += b.y; return a; }
vec2i &operator-=(vec2i &a, const vec2i &b) { a.x -= b.x; a.y -= b.y; return a; }
vec2i operator+(const vec2i &a, const vec2i &b) { return vec2i(a.x + b.x, a.y + b.y); }
bool operator==(const vec2i &a, const vec2i &b) { return a.x == b.x && a.y == b.y; }
ostream &operator<<(ostream &o, const vec2i &v) { return o << v.x << ", " << v.y; }

class Board {
public:
	char b[BOARDSIZE*BOARDSIZE];
	Board() { memset(b, 0, sizeof(b)); }
	Board(const Board &board) { memcpy(b, board.b, sizeof(b)); }
	Board(const char *str) {
		memset(b, 0, sizeof(b));
		memcpy(b, str, sizeof(b));
	}

	Board &operator=(const Board &board) { memcpy(b, board.b, sizeof(b)); return *this; }
	char &operator() (int i, int j) { return b[i + j * BOARDSIZE]; }
	const char &operator() (int i, int j) const { return b[i + j * BOARDSIZE]; }
	char &operator() (const vec2i &v) { return b[v.x + v.y * BOARDSIZE]; }
	const char &operator() (const vec2i &v) const { return b[v.x + v.y * BOARDSIZE]; }

	static bool inside(const vec2i &v) { return v.x >= 0 && v.y >= 0 && v.x < BOARDSIZE && v.y < BOARDSIZE; }


	bool checkAround(const vec2i &v) const {
		static vec2i ofs[4] = {vec2i(-1,0), vec2i(1,0), vec2i(0,-1), vec2i(0,1)};
		for (int l = 0; l < (int)numberof(ofs); l++) {
			vec2i v2 = v + ofs[l];
			if (!inside(v2)) continue;
			if ((*this)(v2) != '.') return true;	//found a touching thing!
		}
		return false;
	}

	//calculate the volume to surface area ratio
	//the greater this value the more compact the board is
	//a sphere has maximal volume to surface area, mind you.  in R^n at least.
	//in manhattan space (which we're dealing with) its a box
	float calcDensity() const {
		int vol = 0, sa = 0;
		for (int j = 0; j < BOARDSIZE; j++) {
			for (int i = 0; i < BOARDSIZE; i++) {
				if ((*this)(i,j) != '.') {
					vol++;
					//test all four directions for an exposed surface
					//this relies on left-to-right condition evaluation
					if (i == 0 || (*this)(i-1,j) == '.') sa++;
					if (j == 0 || (*this)(i,j-1) == '.') sa++;
					if (i == BOARDSIZE-1 || (*this)(i+1,j) == '.') sa++;
					if (j == BOARDSIZE+1 || (*this)(i,j+1) == '.') sa++;
				}
			}
		}
		assert(vol);	//we better have played something
		return (float)vol / (float)sa;
	}
};
ostream &operator<<(ostream &o, const Board &b) {
	for (int j = 0; j < BOARDSIZE; j++) {
		for (int i = 0; i < BOARDSIZE; i++) {
			o << b(i,j);
		}
		o << endl;
	}
	return o;
}


#ifdef WIN32
class Timer {
public:
	void init() {}
	double getTime() { return (double)clock() / (double)CLOCKS_PER_SEC; }
};
#else
class Timer {
	static struct timezone tz;
	struct timeval firstTime;
public:
	void init() { gettimeofday(&firstTime, &tz); }

	double getTime(void) {
		struct timeval tv;
		gettimeofday(&tv, &tz);
		return
			((double)tv.tv_sec + (double)tv.tv_usec / 1000000.0) -
			((double)firstTime.tv_sec + (double)firstTime.tv_usec / 1000000.0);
	}
};
struct timezone Timer::tz = {0,0};
#endif
Timer timer;

class Bonus {
protected:
	Board bonus;
public:
	//map these as follows: 1=dl 2=tl, 3=dw, 4=tw
	Bonus() : bonus(
		"4..1...4...1..4"
		".3...2...2...3."
		"..3...1.1...3.."
		"1..3...1...3..1"
		"....3.....3...."
		".2...2...2...2."
		"..1...1.1...1.."
		"4..1...3...1..4"
		"..1...1.1...1.."
		".2...2...2...2."
		"....3.....3...."
		"1..3...1...3..1"
		"..3...1.1...3.."
		".3...2...2...3."
		"4..1...4...1..4") {}

public:
	int letterMul(const vec2i &v) {
		assert(Board::inside(v));
		char b = bonus(v);
		if (b == '2') return 3;
		if (b == '1') return 2;
		return 1;
	}

	int wordMul(const vec2i &v) {
		assert(Board::inside(v));
		char b = bonus(v);
		if (b == '4') return 3;
		if (b == '3') return 2;
		return 1;
	}
};
Bonus bonus;

//0-25, add 'a' to find the right score (stay in lcase huh)
int letterScoreTable[26] = {
1,3,3,2,1,4,2,4,1,8,
5,1,3,1,1,3,10,1,1,1,
1,4,4,8,4,10
};

int letterScore(char l) {
	if (l >= 'a' && l <= 'z') return letterScoreTable[l - 'a'];
	return 0;
}

template<int CAPACITY>
class IterVec {
public:
	int vec[CAPACITY];
	int size;

	IterVec() : size(0) { memset(vec, 0, sizeof(vec)); }
	IterVec(const IterVec<CAPACITY> &i) { memcpy(vec, i.vec, sizeof(vec)); size = i.size; }

	void reset(int _size) {
		size = _size;
		assert(size);
		memset(vec, 0, sizeof(int) * size);
		vec[0] = -1;	//for the first test/inc
	}

	bool iter(int sup) {
		int k = 0;
		for(; k < size; k++) {
			vec[k]++;
			if (vec[k] < sup) break;
			vec[k] = 0;
		}
		return k != size;
	}

	bool repeats() {
		for (int i = 1; i < size; i++) {
			for (int j = 0; j < i; j++) {
				if (vec[i] == vec[j]) return true;
			}
		}
		return false;
	}
};
template<int CAPACITY> 
ostream &operator << (ostream &o, IterVec<CAPACITY> &iter) {
	o << "[";
	for (int i = 0; i < iter.size; i++) o << iter.vec[i] << " ";
	return o << "]";
}

char *rtrim(char *line) {
	size_t len = strlen(line);
	while(len) {
		if (line[len-1] == 10) line[len-1] = 0;
		else if (line[len-1] == 13) line[len-1] = 0;
		else break;
	}
	return line;
}

// case-independent (ci) string less_than
// returns true if s1 < s2
struct ci_less : binary_function<string, string, bool> {

	// case-independent (ci) compare_less binary function
	struct nocase_compare : public binary_function<unsigned char,unsigned char,bool> {
		bool operator() (const unsigned char& c1, const unsigned char& c2) const {
			return tolower (c1) < tolower (c2); }
	};

	bool operator() (const string & s1, const string & s2) const {
		return lexicographical_compare 
			(s1.begin (), s1.end (),   // source range
			s2.begin (), s2.end (),   // dest range
			nocase_compare ());  // comparison
	}
}; // end of ci_less

//sort them by their length
class Dictionary {
public:
	bool useDict;
	set<string, ci_less> dict;		//if only it could map to n-nested maps ... so its all log(n) access

	Dictionary() : useDict(false) {}

	bool init(const string &filename) {
		ifstream file(filename.c_str());
		if (!file.is_open()) {
			cerr << "failed to find dictionary file " << filename << endl;
			return false;
		}
		while (!file.eof()) {
			char line[256];
			file.getline(line, sizeof(line));
			rtrim(line);
			if (line[0]) dict.insert(line);	//so lookups will be faster
		}
		file.close();
		useDict = true;
		return true;
	}

	//returns true if its in the dictionary
	//modifies the word accordingly if it has any spaces in it
	bool hasWord(char *w) { return !useDict || dict.find(w) != dict.end(); }
};
Dictionary dict;

//make this combo once we evaluate a score
class BoardAndScore {
public:
	Board board;
	int score;

	//debug stuff: how we placed the move
	vec2i pos;
	int dir, len;
	IterVec<MAXHANDSIZE> iter;
	int ident;	//the identity of the thread who calculated this move
	float density;

	BoardAndScore(const Board &_board, int _score, 
		const vec2i &_pos, int _dir, int _len, IterVec<MAXHANDSIZE> &_iter, int _ident,
		float _density
	) : board(_board), score(_score),
		pos(_pos), dir(_dir), len(_len), iter(_iter), ident(_ident),
		density(_density)
	{}

	void report() {
		cout << "pos " << pos << " dir " << dir << " len " << len << " iter " << iter << " ident " << ident << " density " << density << endl;
		cout << "score is " << score << endl;
		cout << board;
		cout << endl;
	}
};

typedef binary_function<BoardAndScore*, BoardAndScore*, bool> ci_bs;

struct ci_score : ci_bs {
	bool operator() (const BoardAndScore *s1, const BoardAndScore *s2) const {
		return s1->score < s2->score;
	}
};

struct ci_density : ci_bs {
	bool operator() (const BoardAndScore *s1, const BoardAndScore *s2) const {
		return s1->density < s2->density;
	}
};


typedef vector<BoardAndScore*> bsvec_t;

class ThreadSafe {
protected:
#if NUMBER_OF_THREADS == 1
	void lock() {}
	void unlock() {}
#else
#ifdef WIN32
	CRITICAL_SECTION cs;	//critical sections are mutexes for single processes
	void lock() { EnterCriticalSection(&cs); }
	void unlock () { LeaveCriticalSection(&cs); }
public:
	ThreadSafe() { InitializeCriticalSection(&cs); }
	~ThreadSafe() { DeleteCriticalSection(&cs); }
#else
	pthread_mutex_t mutex;
	void lock() { pthread_mutex_lock(&mutex); }
	void unlock() { pthread_mutex_unlock(&mutex); }
public:
	ThreadSafe() {
		int result = pthread_mutex_init(&mutex, NULL);
		assert(!result);
	}
	~ThreadSafe() { pthread_mutex_destroy(&mutex); }
#endif
#endif
};

//then use a custom comparator class with a set to sort our boards by score
//then we can list them all, worse-to-best
class ScoreSet : public ThreadSafe {
protected:
	//storing these as a set is keen and all, but it will overwrite boards of teh same score
	//so il'l just sort it at the end
	bsvec_t bs;

public:
	bool sortByDensity;

	ScoreSet() : sortByDensity(false) {}
	
	~ScoreSet() {
		for (size_t i = 0; i < bs.size(); i++) {
			assert(bs[i]);
			delete bs[i];
		}
	}

	void add(const Board &b, int s,
		const vec2i &pos, int dir, int len, IterVec<MAXHANDSIZE> &iter, int ident,
		float density
	) {
		lock();
		BoardAndScore *bas = new BoardAndScore(b,s, pos,dir,len,iter,ident,density);
		assert(bas);
		bs.push_back(bas);
		unlock();
	}

	void report(bool showAll) {
		lock();
		if (showAll) {	//showall means sort then display each in order
			if (sortByDensity) {
				sort<bsvec_t::iterator, ci_density>(bs.begin(), bs.end(), ci_density());
			} else {
				sort<bsvec_t::iterator, ci_score>(bs.begin(), bs.end(), ci_score());
			}
			for (bsvec_t::iterator i = bs.begin(); i != bs.end(); i++) {
				(*i)->report();
			}
		} else {		//otherwise just linearly search through
			BoardAndScore *best = NULL;
			for (bsvec_t::iterator i = bs.begin(); i != bs.end(); i++) {
				if (!best) {
					best = *i;
				} else {
					if (sortByDensity) {
						if (ci_density()(best, *i)) {
							best = *i;
						}
					} else {
						if (ci_score()(best, *i)) {
							best = *i;
						}
					}
				}
			}
			if (best) best->report();
		}
		unlock();
	}
};

//different movement vectors for each direction
vec2i dirmove[2] = {
	vec2i(1,0),
	vec2i(0,1)
};

string hand;	//arbitrary sized, but must be 1<=size<=MAXHANDSIZE
Board board;		//is there a way to make this a read-only board?  so it only gets info on ctor?
bool isStart = true;	//true by default.  false set on prog startup
double startTime = 0.0;	//set on startup
ScoreSet scoreset;

//pass it a vector of played places
//use the dif between board and testboard to figure out whats what
//returns -1 on fail
int getScore(const vec2i *vs, int vsize, const Board &testboard) {
	int totalscore = 0;

	//for each dir
	for (int dir = 0; dir < 2; dir++) {
		Board playboard;	//for each dir clear the board for keeping track of moves

		//for each location played
		for (int i = 0; i < vsize; i++) {
			vec2i v = vs[i];

			//cycle as far back (left/up) as you can
			while(testboard(v) != '.' && Board::inside(v)) v -= dirmove[dir];
			v += dirmove[dir];	//step back in the board

			//see if that word has been recorded (o1 if we make a separate board for this)
			if (playboard(v)) continue;
			//if it has then skip it

			//otherwise push it on the queue, count up its score, and add that to the total
			playboard(v) = 1;

			//now for the score ...
			//actually we should word test too while we're here
			char word[BOARDSIZE+1];
			int wordsize = 0;
			memset(word,0,sizeof(word));

			int wordscore = 0;
			int wordmul = 1;

			while(Board::inside(v) && testboard(v) != '.') {
				char p = testboard(v);
				word[wordsize] = p;
				wordsize++;
				int ls = letterScore(p);		//get the letter score
				if (board(v) == '.') {			//if we're playing this piece
					ls *= bonus.letterMul(v); 	//apply any letter score scales
					wordmul *= bonus.wordMul(v);	//and accum all word score scales
				}
				wordscore += ls;			//then add it to the word score
				v += dirmove[dir];
			}

			if (wordsize == 1) continue;

			if (!dict.hasWord(word)) return -1;

			wordscore *= wordmul;
			totalscore += wordscore;
		}
	}
	//now here's an interesting question.
	//do you get 50 points for playing all seven tiles, or for emptying your hand?
	//more specifically, on the last play if we empty our hands do we get the bonus?
	if (vsize == MAXHANDSIZE) totalscore += 50;
	return totalscore;
}

//board - the base board
//testboard - the board to test the score of
//bestboard - what to overwrite if its better
void testBest(const vec2i *placeloc, int placelocsize, const Board &testboard,
	//debugging.  why so many repeated solutions?
	const vec2i &pos, int dir, int len, IterVec<MAXHANDSIZE> &iter, int ident
) {
	//now for scores
	//for each placeloc
	int score = getScore(placeloc, placelocsize, testboard);
	if (score != -1) {
		scoreset.add(testboard,score,
			pos,dir,len,iter,ident, testboard.calcDensity());
	}
}

void testAllMoves(const vec2i &v, int dir, int l, int ident) {

	//now no matter what our options, if this location starts on a tile then skip it
	//leave it for another length and another position to test
	if (board(v) != '.') return;

	//ok start at i,j and place along 'dir' as long as we can go
	//if we hit a square then skip it
	//if we reach the end of the world then bail out
	//(maybe i should check all lengths as i go or something? less repetition)
	//as we play, check the original board for touching tiles
	//if - by the end - none are found then bail out

	//for all permutations
	IterVec<MAXHANDSIZE> iter;
	for (iter.reset(l); iter.iter((int)hand.length()); ) {
		//make sure there are no repeats
		if (iter.repeats()) continue;
		//iter[0..l-1] holds our permutation
		//now for the tiles given, simulate a placement with our word

		//and cross check it to the dictionary.  that'll be slow.
		Board testboard(board);

		vec2i placeloc[MAXHANDSIZE];
		int placelocsize = 0;

		vec2i v2(v);
		bool touches = false;
		for (int k = 0; k < l; k++) {
			//find an empty location
			while(testboard(v2) != '.') {
				v2 += dirmove[dir];
				if (!Board::inside(v2)) break;
			}
			if (!Board::inside(v2)) break;

			//play there
			char letter = hand[iter.vec[k]];
			testboard(v2) = letter;

			placeloc[placelocsize] = v2;
			placelocsize++;
			//and do the step? 
			//nah, that'll be covered by the loop when it sees the newly placed tile there

			//make sure we're touching prev played words
			//only need to find if its touching once
			//also, dont bother hunt if we're starting.  we operate on dif constraints then.
			if (!touches && !isStart) touches = board.checkAround(v2); 
		}
		if (isStart) {
			if (testboard(7,7) == '.') continue;	//starting boards need to go thru H8
		} else {
			if (!touches) continue;	//wasnt touching anything
		}
		if (!Board::inside(v2)) continue;	//went over the edge!

		//here, for each piece played
		//if any are blanks then iterate through all substitutions
		//make the substitution and then get/test the score
		{
			vec2i spacelocs[MAXHANDSIZE];
			int numspaces = 0;
			for (int i = 0; i < placelocsize; i++) {
				if (testboard(placeloc[i]) == ' ') {
					spacelocs[numspaces] = placeloc[i];
					numspaces++;
				}
			}

			if (!numspaces) {	//if none then just test the original board
				testBest(placeloc, placelocsize, testboard, v, dir, l, iter, ident);
			} else {		//else test all combinations

				//make a board for substituting our blanks with hack letters
				Board blankboard(testboard);

				//now for each space (n-nested loops) 
				IterVec<MAXHANDSIZE> iter;
				for (iter.reset(numspaces); iter.iter(26); ) {
					//fill in the spaces
					for (int i = 0; i < numspaces; i++) {
						//use upper case to denote that it is a blank
						blankboard(spacelocs[i]) = iter.vec[i] + 'A';    
					} 
					testBest(placeloc, placelocsize, blankboard, v, dir, l, iter, ident);
				}
			}
		}
		
		//give us an update every second
		static int lasts = -1;
		int s = (int)(timer.getTime() - startTime);
		if (s != lasts) {
			lasts = s;
			//cerr for info output not for the boards
			cerr << "thread " << ident;
			cerr << " dir = " << dir;
			cerr << " pos = " << v;
			cerr << " iter = " << iter;
			cerr << endl;
		}
	}
}

class BoardMove {
public:
	int l;		//length of the hand, from 1 to hand.length()
	int dir;	//direction we're playing
	vec2i v;	//position in the board, from 0 to BOARDSIZE-l

	BoardMove() : l(1), dir(0), v(0,0) {}
	BoardMove(const BoardMove &s) : l(s.l), dir(s.dir), v(s.v) {}

	BoardMove &operator=(const BoardMove &s) { l = s.l; dir = s.dir; v = s.v; return *this; }

	//increment the iterator
	void inc() {
		if (done()) return;

		v.y++;
		int maxy = BOARDSIZE - dirmove[dir].y * (l-1);
		if (v.y < maxy) return;
		v.y = 0;

		v.x++;
		int maxx = BOARDSIZE - dirmove[dir].x * (l-1);
		if (v.x < maxx) return;
		v.x = 0;

		dir++;
		int maxdir = 2;	//maxdir is usually 2: horizontal vs. vertical
		if (l == 1) maxdir = 1;	//but if we're only playing one tile then direction doesn't matter
		if (dir < maxdir) return;
		dir = 0;

		l++;
		if (l <= (int)hand.length()) return;
		//and l >= hand.length() is our done indicator!
	}

	bool done() { return l > (int)hand.length(); }
};

//has a lot in common with ScoreSet.  consider superclassing.
class SearchState : ThreadSafe {
protected:
	BoardMove ss;
public:

	BoardMove cur() {
		lock();
		BoardMove r(ss);
		unlock();
		return r;
	}

	BoardMove next() {
		lock();
		BoardMove r(ss);
		unlock();
		ss.inc();
		return r;
	}
	bool done() { return ss.done(); }
};
SearchState searchstate;

class ThreadInfo {
public:
	int ident;
	ThreadInfo() : ident(0) {}
};

//i think they might just double.  so long as void* and lpvoid* work the same.  meh.
#ifdef WIN32
DWORD WINAPI searchThread(LPVOID _ptr)
#else
void *searchThread(void *_ptr)
#endif
{
	ThreadInfo *info = (ThreadInfo*)_ptr;
	assert(info);
	for(;;) {
		BoardMove ss = searchstate.next();
		if (ss.done()) break;
		testAllMoves(ss.v, ss.dir, ss.l, info->ident);
	}
	return 0;
}

int main(int argc, char **argv) {
	if (argc < 3) {
		cerr << "usage: scrabble <boardfile> <hand>" << endl;
		return 1;
	}
	const char *filename = argv[1];

	//init the hand
	{
		const char *handstr = argv[2];
		size_t len = strlen(handstr);
		if (len < 1 || len > MAXHANDSIZE) {
			cerr << "you only provided a hand of " << len << " letters" << endl;
			return 1;
		}
		//hands provided will be lowercase'd.  other chars will be assumed to be spaces
		char handarr[MAXHANDSIZE+1];
		memset(handarr, 0, sizeof(handarr));
		for (size_t i = 0; i < strlen(handstr); i++) {
			if (handstr[i] >= 'a' && handstr[i] <= 'z') handarr[i] = handstr[i];
			else if (handstr[i] >= 'A' && handstr[i] <= 'Z') handarr[i] = handstr[i] - 'A' + 'a';
			else handarr[i] = ' ';	//space is our blanks.
		}
		//copy the converted buffer over
		hand = handarr;
	}

	//start off our board
	string boardstr;
	{
		ifstream file(filename);
		if (!file.is_open()) {
			cerr << "failed to load the file " << filename << endl;
			return 1;
		}
		ostringstream s;
		char line[256];
		while (file.getline(line, sizeof(line))) s << rtrim(line);
		boardstr = s.str();
		file.close();
	}
	board = Board(boardstr.c_str());
	isStart = true;	//just in case
	for (int i = 0; i < BOARDSIZE; i++) {
		for (int j = 0; j < BOARDSIZE; j++) {
			//boards are going to be funny
			//letters 'a'-'z' are letters
			//'A'-'Z' are going be be assumed to be blanks switched to that kind of letter
			//correct cases
			if (board(i,j) >= 'A' && board(i,j) <= 'Z') {} //blanks 
			else if (board(i,j) >= 'a' && board(i,j) <= 'z') {} //nothing
			else if (board(i,j) == '.') {} //nothing / preserve
			else {
				cerr << "found an unknown char in the board: " << board(i,j) << " / (" << (int)board(i,j) << ")" << endl;
				return 1;
			}
			if (board(i,j) != '.') {
				isStart = false;
				break;
			}
		}
	}

	bool showAll = false;
	bool useDict = true;
	string dictFilename = "dictionary.txt";

	//extra params
	for (int i = 3; i < argc; i++) {
		//zero-param
		if (!strcasecmp(argv[i], "density")) { scoreset.sortByDensity = true; }
		else if (!strcasecmp(argv[i], "nodict")) { useDict = false; }
		else if (!strcasecmp(argv[i], "showall")) { showAll = true; }
		//one-param
		else if (i < argc-1) {
			if (!strcasecmp(argv[i], "dict")) { dictFilename = argv[i+1]; i++; }
		}
	}

	//load the dictionary
	if (useDict && !dict.init(dictFilename)) return 1;
	cout << "using dictionary? " << (dict.useDict ? dictFilename.c_str() : "/no/") << endl;

	//start your engines
	timer.init();

	startTime = timer.getTime();

#if NUMBER_OF_THREADS == 1
	ThreadInfo threadArgs;
	threadArgs.ident = 0;
	searchThread(&threadArgs);
#else
	//start off the threads!
	ThreadInfo threadArgs[NUMBER_OF_THREADS];	//its either this or find how to poll thread done's from linux

#ifdef WIN32
	HANDLE threads[NUMBER_OF_THREADS];
#else
	pthread_t threads[NUMBER_OF_THREADS];
#endif
	for (int i = 0; i < NUMBER_OF_THREADS; i++) {
		threadArgs[i].ident = i;
#ifdef WIN32
		unsigned long threadID;	//what is this used for? thread identifier?
		threads[i] = CreateThread(NULL, 0, searchThread, threadArgs+i, 0, &threadID);
		assert(threads[i]);
#else
		pthread_create(threads+i, NULL, searchThread, threadArgs+i);
#endif
	}

	for (int i = 0; i < NUMBER_OF_THREADS; i++) {
#ifdef WIN32
		DWORD result = WaitForSingleObject(threads[i], INFINITE);
		assert(result != WAIT_FAILED);
#else
		pthread_join(threads[i], NULL);
#endif
	}
#endif

	//calc time before dumping report
	double timeTaken = timer.getTime() - startTime;

	scoreset.report(showAll);
	cerr << "took " << timeTaken << " seconds to cycle through all those tiles and permutations ..." << endl;

	return 0;
}

