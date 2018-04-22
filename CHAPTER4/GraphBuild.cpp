#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

#define N 10000;

class Edge {

public:

	int st;
	int end;
	int weight;

	Edge() {}

	Edge(int st, int end, int weight) {
		this->st = st;
		this->end = end;
		this->weight = weight;
	}

	~Edge() {}

};

class MinHeap {

public:

	Edge * HeapArr;
	int curCount;

	MinHeap();
	~MinHeap() { delete[] HeapArr; }

	void Insert(Edge eg);
	void SiftDown(int inx);
	void SiftUp();
	Edge RemoveMin();
};

MinHeap::MinHeap() {

	HeapArr = new Edge[36];
	curCount = 0;

}

void MinHeap::SiftDown(int inx) {

	
	Edge edgeTmp = HeapArr[inx];
	int j = 2 * inx + 1;

	while (j < curCount) {
		if (j< curCount - 1 && HeapArr[j].weight> HeapArr[j + 1].weight)
			j++;

		if (HeapArr[j].weight < edgeTmp.weight) {
			HeapArr[inx] = HeapArr[j];
			inx = j;
			j = 2 * j + 1;
		}
		else
			break;
	}
	HeapArr[inx] = edgeTmp;
}

void MinHeap::Insert(Edge eg) {

	HeapArr[curCount] = eg;
	SiftUp();
	curCount++;
}

void MinHeap::SiftUp() {

	int inx = this->curCount;
	int tmp = -1;
	Edge edgeTmp = HeapArr[inx];

	while (inx != 0) {

		tmp = inx / 2 - 1 + (inx % 2);
		if (HeapArr[inx].weight < HeapArr[tmp].weight) {
			HeapArr[inx] = HeapArr[tmp];
			inx = tmp;
		}

		else
			break;
	}

	HeapArr[inx] = edgeTmp;
}

Edge MinHeap::RemoveMin() {

	Edge min = HeapArr[0];
	HeapArr[0] = HeapArr[curCount - 1];
	curCount--;

	for (int i = curCount / 2 - 1; i >= 0; i--)
		SiftDown(i);

	return min;

}


class ArrStack {

public:

	int top;
	int maxSize;
	int *arr;

	ArrStack(int maxSize) {
		this->maxSize = maxSize;
		top = -1;
		arr = new int[this->maxSize];
	}

	~ArrStack() { delete[] arr; }

	bool Push(int x) {
		if (top == maxSize - 1) {
			cout << "FULL Stack!" << endl;
			return false;
		}

		arr[++top] = x;
		return true;
	}

	int Pop() {
		if (top == -1) {
			cout << "EMPTY Stack!" << endl;
			return -16;
		}

		return arr[top--];
	}

	int Top() {
		if (top == -1) {
			cout << "EMPTY Stack!" << endl;
			return -16;
		}

		return arr[top];
	}

	bool isEmpty() { return (top == -1) ? true : false; }

};


class ArrQueue {

public:

	int rear;
	int front;
	int maxSize;
	int *q;

	ArrQueue(int maxSize) {
		front = rear = 0;
		this->maxSize = maxSize;
		q = new int[this->maxSize + 1];
	}

	~ArrQueue() { delete[] q; }

	bool EnQueue(int x) {
		if ((rear + 1) % this->maxSize == front) {
			cout << "FULL Queue!" << endl;
			return false;
		}

		q[rear] = x;
		rear = (rear + 1) % this->maxSize;
		return true;
	}

	int DeQueue() {
		if (isEmpty()) {
			cout << "EMPTY Queue!" << endl;
			return false;
		}

		int res = q[front];
		front = (front + 1) % this->maxSize;
		return res;
	}

	int GetFront() {
		if (isEmpty()) {
			cout << "EMPTY Queue!" << endl;
			return false;
		}

		return q[front];
	}

	bool isEmpty() { return (rear == front) ? true : false; }

};


class UFSets {

public:

	int size;
	int *root;
	int *next;
	int *length;

	UFSets(int size) {
		this->size = size;
		root = new int[size];
		next = new int[size];
		length = new int[size];

		for (int i = 0; i < size; i++) {
			root[i] = next[i] = i;
			length[i] = 1;
		}
	}


	void UFSetsUnion(int v, int u) {
		if (root[u] == root[v])
			return;

		length[v] < length[u] ? Union(v, u) : Union(u, v);
	}

	void Union(int varA, int varB) { // varA with shorter length

		int rt = root[varA];
		length[root[varB]] += length[root[varA]];

		root[rt] = root[varB];
		for (int i = next[rt]; i != rt; i = next[i])
			root[i] = root[varB];

		Swap(varA, varB);

	}

	void Swap(int varA, int varB) {

		int nextA = next[varA];
		int nextB = next[varB];

		next[varA] = nextB;
		next[varB] = nextA;
	}

};


class AdjGraph {

public:

	Edge * EdgeArr;
	int ** AdjMatrix;

	int VertexNum;
	int EdgeNum;

	AdjGraph(int verNum, int opt);
	AdjGraph(int verNum);
	~AdjGraph() {}

	void SetEdge(int st, int end); // no direction
	void SetEdge(int st, int end, int weight); // no direction, record in EdgeArr
	void SetFTEdge(int st, int end, int weight) { this->AdjMatrix[st][end] = weight; } // with direction

	Edge FindFristEdge(int vertex);
	Edge FindFristEdge(int vertex, int back);
	Edge FindNextEdge(Edge oneEdge);

	void BFSTraverse();
	void BFS(int inx, bool mark[]);

	void DFSNoReverse();
	void DFS(int inx, bool mark[]);
	void DFSTraverse();

	Edge * Prim();
	Edge * Kruskal();

	void Dijkstra(int inx);
	void Floyd();

	int DFSToCircle();
	void FindCircle();
	int DFS4Circle(int inx, bool isPassed[], int circlePath[], int prevPath[]);
	Edge* BreakCircle();

};

Edge* AdjGraph::BreakCircle() {

	Edge* resEdge = new Edge[this->VertexNum - 1];

	int stPos = -1, endPos= -1, maxInxSt= -1, maxInxEnd = -1, maxWeight=-1, stMark = -1;
	bool* isPassed = new bool[this->VertexNum];
	int* circlePath = new int[this->VertexNum];
	int* prevPath = new int[this->VertexNum];

	for (int i = 0; i < this->VertexNum; i++) {
		isPassed[i] = false;
		circlePath[i] = -1;
		prevPath[i] = -1;
	}


	for (int i = this->EdgeNum; i > this->VertexNum - 1; i--) {

		stPos = DFS4Circle(0, isPassed, circlePath, prevPath);
		stMark = stPos;
		endPos = circlePath[stPos];
		
		maxWeight = this->AdjMatrix[stPos][endPos];
		maxInxSt = stPos;
		maxInxEnd = endPos;
		

		for (endPos; endPos != stMark;endPos = circlePath[endPos]) {

			if (this->AdjMatrix[endPos][circlePath[endPos]] > maxWeight) {
				maxWeight = this->AdjMatrix[endPos][circlePath[endPos]];
				maxInxSt = endPos;
				maxInxEnd = circlePath[endPos];
			}
		}

		this->AdjMatrix[maxInxSt][maxInxEnd] = this->AdjMatrix[maxInxEnd][maxInxSt] = 0;


		for (int i = 0; i < this->VertexNum; i++) {
			isPassed[i] = false;
			circlePath[i] = -1;
			prevPath[i] = -1;
		}
	}

	for (int i = 0; i < this->VertexNum; i++) {
		for (int j = 0; j < this->VertexNum; j++)
			cout << this->AdjMatrix[i][j] << ' ';
		cout << endl;
	}

	return NULL;
}

int AdjGraph::DFS4Circle(int inx, bool isPassed[], int circlePath[], int prevPath[]) {

	isPassed[inx] = true;
	bool isFound = false;
	int res = -1;

	for (int i = 0; i < this->VertexNum; i++) {
		if (i != prevPath[inx] && this->AdjMatrix[inx][i] != 0) {
			prevPath[i] = inx;
			circlePath[inx] = i;
			isFound = true;
			if (isPassed[i])
				return i;
			res = DFS4Circle(i, isPassed, circlePath, prevPath);
			if (res >= 0)
				return res;
			else
				isFound = false;
		}
	}

	if (!isFound)
		return -1;
}

/*
Edge* AdjGraph::BreakCircle() {

	Edge* resEdge = new Edge[this->VertexNum - 1];
	int stPos = DFSToCircle(), endPos = -1, tmp = -1;
	int stMark = -1, maxWeight = -1, maxInxSt = -1, maxInxEnd = -1;

	bool * visMark = new bool[this->VertexNum];
	for (int i = 0; i < this->VertexNum; i++)
		visMark[i] = false;

	bool isBack = false;


	while (stPos != -1) {
		stMark = stPos;
		endPos = FindFristEdge(stPos).end;
		maxInxSt = stPos;
		maxInxEnd = endPos;
		maxWeight = this->AdjMatrix[maxInxSt][maxInxEnd];

		tmp = stPos;
		stPos = endPos;
		endPos = FindFristEdge(endPos, tmp).end;
		while (true) {
			
			visMark[stPos] = true;

			if (this->AdjMatrix[stPos][endPos] > maxWeight) {
				maxInxSt = stPos;
				maxInxEnd = endPos;
				maxWeight = this->AdjMatrix[maxInxSt][maxInxEnd];
			}

			if (isBack)
				break;

			if (this->AdjMatrix[endPos][stMark] != 0) {
				isBack = true;
				stPos = endPos;
				endPos = stMark;
				
			}

			else {
				stPos = endPos;
				//visMark[stPos] = true;
				for (int i = 0; i < this->VertexNum; i++)
					if (this->AdjMatrix[stPos][i] != 0 && !visMark[i]) {
						endPos = i;
						break;
					}
			}

		}

		isBack = false;
		
		this->AdjMatrix[maxInxSt][maxInxEnd] = this->AdjMatrix[maxInxEnd][maxInxSt] = 0;
		stPos = DFSToCircle();
		tmp = -1;
		for (int i = 0; i < this->VertexNum; i++)
			visMark[i] = false;
	}

	for (int i = 0; i < this->VertexNum; i++) {
		for (int j = 0; j < this->VertexNum; j++)
			cout << this->AdjMatrix[i][j] << ' ';
		cout << endl;
	}

	return NULL;

}
*/
/*
int AdjGraph::DFSToCircle() {

	ArrStack as(36);
	bool *visMark = new bool[this->VertexNum];
	int *prevPath = new int[this->VertexNum];
	bool *isPassed = new bool[this->VertexNum];
	int curPos = 0, prev = -1, resSt = -1;
	for (int i = 0; i < this->VertexNum; i++) {
		visMark[i] = false;
		prevPath[i] = -1;
		isPassed[i] = false;
	}

	bool isFound = false;

	for (int inx = 0; inx < this->VertexNum; inx++) {
		if (!visMark[inx]) {
			as.Push(inx);
			visMark[inx] = true;
		}
		while (!as.isEmpty()) {
			curPos = as.Pop();
			isPassed[curPos] = true;
			for (int i = 0; i < this->VertexNum; i++) {

				if (this->AdjMatrix[curPos][i] != 0 && isPassed[i] && i!= prevPath[curPos]) { //  i != prev && 
					isFound = true;
					resSt = i;
					break;
				}

				if (this->AdjMatrix[curPos][i] != 0 && !visMark[i]) {
					as.Push(i);
					visMark[i] = true;
					prevPath[i] = curPos;
				}

				if (this->AdjMatrix[curPos][i] != 0 && !isPassed[i] && visMark[i])
					prevPath[i] = curPos;

			}
			if (isFound)
				break;
		}
		if (isFound)
			break;
	}
	if (isFound)
		return resSt;
	else
		return -1;
}
*/
void AdjGraph::FindCircle() {

	bool * visMark = new bool[this->VertexNum];
	int * inDegree = new int[this->VertexNum];
	for (int i = 0; i < this->VertexNum; i++) {
		visMark[i] = false;
		inDegree[i] = 0;
	}

	bool isFound = false;
	int pos = 0, Counter = 0;
	

	for (int i = 0; i < this->VertexNum; i++)
		for (int j = 0; j < this->VertexNum; j++)
			if (this->AdjMatrix[j][i] != 0)
				inDegree[i]++;



	for (int i = 0; i < this->VertexNum; i++) {
		for (pos = 0; pos < this->VertexNum; pos++)
			if (inDegree[pos] == 0 && !visMark[pos]) {
				isFound = true;
				visMark[pos] = true;
				Counter++;
				break;
			}

		if (isFound) {
			for (int j = 0; j < this->VertexNum; j++)
				if (this->AdjMatrix[pos][j] != 0 && !visMark[j]) {
					this->AdjMatrix[pos][j] = 0;
					inDegree[j]--;
				}
			isFound = false;
		}

		else
			break;
	}

	if (Counter == this->VertexNum)
		cout << "There is no Circle!";

	else {
		cout << "There is at least one Circle: ";
		int inx = 0, stPos;

		bool* isPassed = new bool[this->VertexNum];
		int* circlePath = new int[this->VertexNum];
		int* prevPath = new int[this->VertexNum];

		for (int i = 0; i < this->VertexNum; i++) {
			isPassed[i] = false;
			circlePath[i] = -1;
			prevPath[i] = -1;
		}

		for (inx; inx < this->VertexNum; inx++)
			if (!visMark[inx])
				break;
		
		stPos = DFS4Circle(inx, isPassed, circlePath, prevPath);
		cout << stPos << ' ';
		for (int endPos = circlePath[stPos]; endPos != stPos; endPos = circlePath[endPos])
			cout << endPos << ' ';
		cout << stPos << ' ';
		/*
		for (inx; inx < this->VertexNum; inx++)
			if (!visMark[inx])
				break;

		stPos = inx;
		cout << stPos << ' ';
		while (true) {
			for (int i = 0; i < this->VertexNum; i++)
				if (this->AdjMatrix[inx][i] != 0) {
					inx = i;
					cout << inx << ' ';
					break;
				}

			if (inx == stPos)
				break;
		}*/
	}
	cout << endl;


}


void AdjGraph::Dijkstra(int inx) {

	int* disArr = new int[this->VertexNum];
	int* prevPath = new int[this->VertexNum];
	bool* visMark = new bool[this->VertexNum];

	for (int i = 0; i < this->VertexNum; i++) {
		disArr[i] = N;
		prevPath[i] = -1;
		visMark[i] = false;
	}

	visMark[inx] = true;
	disArr[inx] = 0;
	prevPath[inx] = 0;
	int minWeight, minInx, curInx = inx;

	for (int i = 0; i < this->VertexNum - 1; i++) {

		minWeight = 10000;
		minInx = -1;

		for (int j = 0; j < this->VertexNum; j++) {
			if ((disArr[curInx] + this->AdjMatrix[curInx][j]) < disArr[j] && !visMark[j]) {   //this->AdjMatrix[curInx][j] != 0 && 
				disArr[j] = disArr[curInx] + this->AdjMatrix[curInx][j];
				prevPath[j] = curInx;
			}
		}

		for (int j = 0; j < this->VertexNum; j++) {
			if (!visMark[j]) {
				if (disArr[j] < minWeight) {
					minWeight = disArr[j];
					minInx = j;
				}
			}
		}

		curInx = minInx;
		visMark[minInx] = true;
	}

	ArrStack aq(36);

	

	for (int i = 0; i < this->VertexNum; i++) {
		
		cout << inx << " to " << i << ": ";
		if (i != inx)
			aq.Push(i);
		for (int j = i; prevPath[j] != inx; j = prevPath[j])
			aq.Push(prevPath[j]);

		aq.Push(inx);
		while (!aq.isEmpty())
			cout << aq.Pop() << ' ';
		cout << endl;
		
	}

}


void AdjGraph::Floyd() {
	 
	int** adjMat, **prevPath;
	adjMat = new int*[this->VertexNum];
	prevPath = new int*[this->VertexNum];
	for (int i = 0; i < this->VertexNum; i++) {
		adjMat[i] = new int[this->VertexNum];
		prevPath[i] = new int[this->VertexNum];
	}

	for (int i = 0; i < this->VertexNum; i++)
		for (int j = 0; j < this->VertexNum; j++) {
			prevPath[i][j] = i;
			adjMat[i][j] = this->AdjMatrix[i][j];
		}

	for (int i = 0; i < this->VertexNum; i++) {

		for(int j = 0;j<this->VertexNum;j++)
			for (int k = 0; k < this->VertexNum; k++) {
				if (j != i && k != i) {
					if (adjMat[j][k] > adjMat[j][i] + adjMat[i][k]) {
						adjMat[j][k] = adjMat[j][i] + adjMat[i][k];
						prevPath[j][k] = i;
					}
				}
			}
	}

	for (int i = 0; i<this->VertexNum; i++)
	{
		for (int j = 0; j<this->VertexNum; j++)
			cout << prevPath[i][j] << ' ';
		cout << endl;
	}

}



Edge* AdjGraph::Kruskal() {

	UFSets ufs(this->VertexNum);
	MinHeap heap;
	Edge* resEdge = new Edge[this->VertexNum - 1];

	for (int i = 0; i < this->VertexNum; i++) {
		for (int j = i + 1; j < VertexNum; j++) {
			if (this->AdjMatrix[i][j] != 0)
				heap.Insert(Edge(i, j, this->AdjMatrix[i][j]));
		}
	}

	Edge tmpEdge;
	int st, end, weight, Counter = 0;

	while (Counter < 5){

		tmpEdge = heap.RemoveMin();
		st = tmpEdge.st;
		end = tmpEdge.end;

		if (ufs.root[st] != ufs.root[end]) {
			resEdge[Counter++] = tmpEdge;
			ufs.UFSetsUnion(st, end);
		}
	}

	return resEdge;
}


Edge* AdjGraph::Prim(){

	Edge * resEdge = new Edge[this->VertexNum];
	
	int* nearestWeight = new int[this->VertexNum];
	int* neighborInx = new int[this->VertexNum];

	for (int i = 0; i < this->VertexNum; i++) {
		neighborInx[i] = 0;
		nearestWeight[i] = N;

	}
	for (int i = 1; i < this->VertexNum; i++)
		if (this->AdjMatrix[i][0] != 0)
			nearestWeight[i] = this->AdjMatrix[i][0];

	neighborInx[0] = -1;

	for (int i = 1; i < this->VertexNum; i++) {
		int min = 5000, min_inx = -1;
		for (int j = 0; j < this->VertexNum; j++)
			if (nearestWeight[j] < min && neighborInx[j] >-1) {
				min = nearestWeight[j];
				min_inx = j;
			}

		if (min_inx >= 0) {
			resEdge[i - 1] = Edge(neighborInx[min_inx], min_inx, min);
			neighborInx[min_inx] = -1;

			for (int j = 0; j < this->VertexNum; j++) {
				if (neighborInx[j] > -1) {
					if (this->AdjMatrix[min_inx][j] != 0 && this->AdjMatrix[min_inx][j] < nearestWeight[j]) {
						nearestWeight[j] = this->AdjMatrix[min_inx][j];
						neighborInx[j] = min_inx;
					}
				}
			}

		}
	}

	return resEdge;

}


void AdjGraph::DFSNoReverse() {

	ArrStack as(36);
	bool * mark = new bool[this->VertexNum];

	for (int i = 0; i < this->VertexNum; i++)
		mark[i] = false;
	int tmp = -1;

	for(int inx=0;inx<this->VertexNum;inx++){

		if (!mark[inx]) {
			as.Push(inx);
			mark[inx] = true;
		}
		while (!as.isEmpty()) {
			tmp = as.Pop();
			cout << tmp + 1 << ' ';
			//mark[tmp] = true;
			for (int i = 0; i < this->VertexNum; i++) {
				if (this->AdjMatrix[tmp][i] != 0 && !mark[i]) {
					as.Push(i);
					mark[i] = true;
				}
			}
		}
	}

}

void AdjGraph::DFSTraverse() {

	bool * mark = new bool[this->VertexNum];

	for (int i = 0; i < this->VertexNum; i++)
		mark[i] = false;

	for (int i = 0; i < this->VertexNum; i++)
		if (!mark[i])
			DFS(i, mark);

}

void AdjGraph::DFS(int inx, bool mark[]) {

	cout << inx + 1 << ' ';
	mark[inx] = true;
	for (int i = 0; i < this->VertexNum; i++) {
		if (this->AdjMatrix[inx][i] != 0 && !mark[i]) {
			DFS(i, mark);
		}
	}

}

void AdjGraph::BFSTraverse() {

	bool * mark = new bool[this->VertexNum];

	for (int i = 0; i < this->VertexNum; i++)
		mark[i] = false;

	for (int i = 0; i < this->VertexNum; i++)
		if (!mark[i])
			BFS(i, mark);

}


void AdjGraph::BFS(int inx, bool mark[]) {

	ArrQueue aq(36);

	if(!mark[inx])
		aq.EnQueue(inx);
	
	int tmp = -1;

	while (!aq.isEmpty()) {
		tmp = aq.DeQueue();
		cout << tmp + 1 << ' ';
		mark[tmp] = true;

		for(int i=0;i<this->VertexNum;i++){
			if (this->AdjMatrix[tmp][i] != 0 && !mark[i]) {
				aq.EnQueue(i);
				mark[i] = true;
			}
		}
	}

}


Edge AdjGraph::FindFristEdge(int vertex) {

	bool isFound = false;
	Edge resEdge;
	resEdge.st = vertex;

	for (int i = 0; i < this->VertexNum; i++) {
		if (this->AdjMatrix[vertex][i] != 0) {
			resEdge.end = i;
			isFound = true;
			break;
		}
	}

	if (!isFound)
		resEdge.end = -1;
	return resEdge;
}

Edge AdjGraph::FindFristEdge(int vertex, int back) {
	bool isFound = false;
	Edge resEdge;
	resEdge.st = vertex;

	for (int i = 0; i < this->VertexNum; i++) {
		if (this->AdjMatrix[vertex][i] != 0 && i != back) {
			resEdge.end = i;
			isFound = true;
			break;
		}
	}


	if (!isFound)
		resEdge.end = -1;
	return resEdge;
}

Edge AdjGraph::FindNextEdge(Edge oneEdge) {
	bool isFound = false;
	Edge resEdge;
	resEdge.st = oneEdge.st;

	for (int i = oneEdge.end; i < this->VertexNum; i++) {
		if (this->AdjMatrix[oneEdge.st][i] != 0) {
			resEdge.end = i;
			isFound = true;
			break;
		}
	}

	if (!isFound)
		resEdge.end = -1;
	return resEdge;
}

void AdjGraph::SetEdge(int st, int end) {

	this->AdjMatrix[st][end] = 1;
	this->AdjMatrix[end][st] = 1;
}

void AdjGraph::SetEdge(int st, int end, int weight) {
	this->AdjMatrix[st][end] = weight;
	this->AdjMatrix[end][st] = weight;
	this->EdgeArr[EdgeNum] = Edge(st, end, weight);
	this->EdgeNum++;

}





AdjGraph::AdjGraph(int verNum, int opt) {

	this->VertexNum = verNum;
	int i, j;

	EdgeArr = new Edge[36];
	EdgeNum = 0;

	this->AdjMatrix = new int*[this->VertexNum];
	for (i = 0; i < verNum; i++)
		AdjMatrix[i] = new int[this->VertexNum];

	for (i = 0; i < this->VertexNum; i++)
		for (j = 0; j < this->VertexNum; j++)
			this->AdjMatrix[i][j] = N;

	for (i = 0; i < this->VertexNum; i++)
		this->AdjMatrix[i][i] = 0;

}

AdjGraph::AdjGraph(int verNum) {

	this->VertexNum = verNum;
	int i, j;

	EdgeArr = new Edge[36];
	EdgeNum = 0;

	this->AdjMatrix = new int*[this->VertexNum];
	for (i = 0; i < verNum; i++)
		AdjMatrix[i] = new int[this->VertexNum];

	for (i = 0; i < this->VertexNum; i++)
		for (j = 0; j < this->VertexNum; j++)
			this->AdjMatrix[i][j] = 0;

}

void main() {
	/*
	AdjGraph gra(8);
	gra.SetEdge(0, 1);
	gra.SetEdge(0, 2);
	gra.SetEdge(1, 3);
	gra.SetEdge(1, 4);
	gra.SetEdge(3, 7);
	gra.SetEdge(2, 5);
	gra.SetEdge(2, 6);
	gra.SetEdge(5, 6);
	gra.SetEdge(7, 4);

	gra.BFSTraverse();
	cout << endl;
	gra.DFSTraverse();
	cout << endl;
	gra.DFSNoReverse();
	*/
	/*
	AdjGraph gra(6);
	gra.SetEdge(0, 1, 6);
	gra.SetEdge(0, 2, 1);
	gra.SetEdge(0, 3, 5);
	gra.SetEdge(1, 2, 5);
	gra.SetEdge(1, 4, 3);
	gra.SetEdge(2, 3, 5);
	gra.SetEdge(2, 4, 6);
	gra.SetEdge(2, 5, 4);
	gra.SetEdge(3, 5, 2);
	gra.SetEdge(4, 5, 6);
	Edge * EdgeShow;
	EdgeShow = gra.Prim();
	for (int i = 0; i<5; i++)
		cout << EdgeShow[i].st + 1 << ' ' << EdgeShow[i].end + 1 << ' ' << EdgeShow[i].weight << endl;
	cout << endl;
	EdgeShow = gra.Kruskal();

	for (int i = 0; i<5; i++)
		cout << EdgeShow[i].st + 1 << ' ' << EdgeShow[i].end + 1 << ' ' << EdgeShow[i].weight << endl;

	system("pause");
	*/
	/*
	AdjGraph gra(6,-1);


	gra.SetFTEdge(0, 1, 12);
	gra.SetFTEdge(0, 2, 10);
	gra.SetFTEdge(0, 4, 30);
	gra.SetFTEdge(0, 5, 100);
	gra.SetFTEdge(1, 2, 5);
	gra.SetFTEdge(2, 3, 50);
	gra.SetFTEdge(3, 5, 10);
	gra.SetFTEdge(4, 3, 20);
	gra.SetFTEdge(4, 5, 60);
	gra.Dijkstra(0);
	*/
	/*
	AdjGraph gra(5,-1);
	gra.SetFTEdge(0, 1, 1);
	gra.SetFTEdge(1, 3, 2);
	gra.SetFTEdge(0, 3, 10);
	gra.SetFTEdge(0, 4, 100);
	gra.SetFTEdge(3, 4, 100);
	gra.SetFTEdge(3, 2, 1);
	gra.SetFTEdge(2, 4, 2);
	gra.Floyd();
	*/

	
	AdjGraph gra(5);
	gra.SetFTEdge(0, 1, 1);
	gra.SetFTEdge(0, 2, 1);
	gra.SetFTEdge(2, 1, 1);
	gra.SetFTEdge(1, 3, 1);
	gra.SetFTEdge(3, 2, 1);
	gra.SetFTEdge(1, 4, 1);
	
	gra.FindCircle();
	system("pause");
	/*
	AdjGraph gra(8);
	gra.SetEdge(0, 1);
	gra.SetEdge(0, 2);
	gra.SetEdge(1, 3);
	//gra.SetEdge(1, 4);
	gra.SetEdge(3, 7);
	gra.SetEdge(2, 5);
	gra.SetEdge(2, 6);
	gra.SetEdge(5, 6);
	gra.SetEdge(7, 4);

	bool* isPassed = new bool[gra.VertexNum];
	int* prevPath = new int[gra.VertexNum];
	int* circlePath = new int[gra.VertexNum];
	for (int i = 0; i < gra.VertexNum; i++) {
		isPassed[i] = false;
		prevPath[i] = -1;
		circlePath[i] = -1;

	}

	int res = gra.DFS4Circle(0, isPassed, circlePath, prevPath);

	cout<<res<<endl;
	for (int i = 0; i < gra.VertexNum; i++) {
		cout << circlePath[i] << ' ';
	}
	cout<<endl;

	//gra.DFSNoReverse();*/

	
	/*
	AdjGraph gra(6);
	gra.SetEdge(0, 1, 6);
	gra.SetEdge(0, 2, 1);
	gra.SetEdge(0, 3, 5);
	gra.SetEdge(1, 2, 5);
	gra.SetEdge(1, 4, 3);
	gra.SetEdge(2, 3, 5);
	gra.SetEdge(2, 4, 6);
	gra.SetEdge(2, 5, 4);
	gra.SetEdge(3, 5, 2);
	gra.SetEdge(4, 5, 6);
	gra.BreakCircle();
	system("pause");
	*/
}