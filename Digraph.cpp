// Digraph.cpp -- directed graph class
// c. 2017 T. O'Neil

#include "Digraph.hpp"
#include <vector>
#include <map>

const int infinity = 99999999;

// there are two implementations, one with min Heap, and the other just by iterating through the whole graph,
// to use the min Heap, set useHeap to true
const bool useHeap = true;

// Node class used in min Heap
struct HNode
{
    // set the default value to "*" and change it when you have an actual value
    // this way you dont have to worry about setting it for combined nodes
    int index = 0;
    int weight = 0;

    // create a new node from a string and its value
    // you must implement this function
    HNode(const int &s, const int &t)
        : index(s), weight(t) {}

    // overloading comparison operators
    // when comparing two pointers to HNode, make sure to de-reference them first
    // i.e. if (*node_pointer1 < *node_pointer2) ...
    bool operator<(const HNode &h) { return weight < h.weight; }
    bool operator<=(const HNode &h) { return weight <= h.weight; }

    // add the friend keyword to allow access to private members
    friend std::ostream &operator<<(std::ostream &os, const HNode *obj)
    {
        if (obj == nullptr)
        {
            os << "nullptr" << std::endl;
            return os;
        }
        os << "Index: " << obj->index << " Weight: " << obj->weight << std::endl;
        return os;
    }
};

// anonymous namespace to hide helper functions, so that they're not confused with the STL functions
namespace
{

    // get the minimum node of two indices
    int min(std::vector<HNode *> &v, int indexA, int indexB)
    {
        if (v[indexA] < v[indexB])
        {
            return indexA;
        }
        else
        {
            return indexB;
        }
    }

    // swap two indices in the tree
    void swap(std::vector<HNode *> &tree, int indexA, int indexB)
    {
        auto temp = tree[indexA];
        tree[indexA] = tree[indexB];
        tree[indexB] = temp;
    }
}

// min Heap class used to implement a priority queue
class Heap
{
protected:
    // holds the tree
    std::vector<HNode *> tree;

    // index of the last element added
    int position = 0;
    int count = 0; // size of tree I think

public:
    Heap() { tree.resize(1); }

    // add an element to the tree
    void enqueue(HNode *node)
    {
        // the first element is always 0 in the heap vector
        // and the true root is at index 1
        tree.push_back(node);
        ++position;
        ++count;

        // we need to heapify the tree from the root.
        fix_up(position);
    }
    // fix the heap from a specific index up
    void fix_up(const int &index)
    {
        // if the heap is only one element, then no need to heapify because it's a heap by definiton.
        // if index is 0 or lower, then it is out of bounds.
        if (index <= 1)
        {
            return;
        }

        // start at the current index, and go up the tree by finding the parent and swaping until we reach the root.
        for (int currentIndex = index; currentIndex > 0; --currentIndex)
        {
            int parentIndex = currentIndex / 2;

            // if the parent is smaller than the current node, then swap them.
            while (parentIndex > 1 && tree[currentIndex]->weight < tree[parentIndex]->weight)
            {
                swap(tree, currentIndex, parentIndex);
                parentIndex /= 2;
            }
        }
    }

    // remove the smallest element
    HNode *dequeue()
    {
        // if tree has 1 or no elements, then it's empty because we don't count the element at index 0.
        if (tree.size() <= 1)
        {
            return nullptr;
        }

        // take first element and replace it with last element to maintain the heap property.
        HNode *minNode = tree[1];
        tree[1] = tree[position];

        tree.pop_back(); // Remove the last element.
        --position;
        --count;

        // we just put the last element in the root, now we need to heapify it from root.
        fix_down(1);

        return minNode;
    }
    // fix the tree after replacing the smallest element
    void fix_down(const int &index)
    {

        // if the heap has only one node or none (since the first node is on index 1)
        if (tree.size() <= 2)
        {
            return;
        }

        int parent = index;
        while (parent < tree.size())
        {
            int left = parent * 2;
            int right = left + 1;
            // if left index exists, then check if its smaller than parent. If right index exists, then check if its smaller than parent.
            // if either are smaller than parent, then we will swap the smaller with the parent
            if ((left < tree.size() && tree[left] < tree[parent]) || (right < tree.size() && tree[right] < tree[parent]))
            {
                int smaller = min(tree, left, right);
                swap(tree, parent, smaller);
                parent = smaller;
            }
            else
            {
                break;
            }
        }
    }
    HNode top()
    {
        return *tree[1];
    }
    // add the friend keyword to allow access to private members
    friend std::ostream &operator<<(std::ostream &os, const Heap &obj)
    {
        for (int i = 0; i < obj.tree.size(); i++)
        {
            os << obj.tree[i];
        }
        return os;
    }

    bool empty()
    {
        return tree.size() <= 1;
    }

    ~Heap()
    {
        for(HNode* node : tree)
        {
            delete node;
        }
        tree.clear();
        position = 0;
        count = 0;
    }
};

// returns the number of vertices
unsigned int Digraph::noVertices()
{
    return numberOfVertices;
}

// returns the number of edges
unsigned int Digraph::noEdges()
{
    return numberOfEdges;
}

// reset all edges to infinity
void Digraph::resetEdges()
{
    numberOfEdges = 0;
    // iterate through all of the distances and set them to infinity
    for (int i = 0; i < numberOfVertices; i++)
    {
        for (int j = 0; j < numberOfVertices; j++)
        {
            distMatrix[i][j] = infinity;
        }
    }
}

// add an edge by setting the distance
void Digraph::addEdge(int source, int dest, int wt)
{
    ++numberOfEdges;
    distMatrix[source][dest] = wt;
}

// delete an edge by setting the distance to infinity
void Digraph::delEdge(int source, int dest)
{
    distMatrix[source][dest] = infinity;
    distMatrix[dest][source] = infinity;
    --numberOfEdges;
}

// returns 0(false) if the distance is infinity, otherwise returns 1(true)
int Digraph::isEdge(int source, int dest)
{
    return distMatrix[source][dest] != infinity;
}

// returns the distance from source to destination using dijkstra's algorithm
int Digraph::dijkstra(int source, int dest)
{
    if (!useHeap)
    {
        // vector contains the distance from vertex[source] to all other vertices
        std::vector<int> weight(numberOfVertices);
        // Initialize all distances to infinity
        for (int i = 0; i < numberOfVertices; i++)
        {
            weight[i] = infinity;
        }

        // The distance of the vertex that we're starting at is 0
        weight[source] = 0;

        // go through the whole graph
        for (int i = 0; i < numberOfVertices; i++)
        {
            // Find closest vertex to source that has not been visited
            int closestIndex = 0;
            for (int j = 0; j < numberOfVertices; j++)
            {
                if (vertex[j]->getStatus() != VISITED)
                {
                    if (weight[j] < weight[closestIndex])
                    {
                        closestIndex = j;
                    }
                }
            }

            // Finalize closest vertex with the group of visited vertices
            vertex[closestIndex]->setStatus(VISITED);

            // Go through again and update distance if there is a better route to take, and ignore the bad routes.
            for (int i = 0; i < numberOfVertices; i++)
            {
                // If there is a distance smaller than current distance[i]
                if (distMatrix[closestIndex][i] != infinity && weight[i] > weight[closestIndex] + distMatrix[closestIndex][i])
                {
                    weight[i] = weight[closestIndex] + distMatrix[closestIndex][i];
                }
            }
        }

        // Reset all vertices to not visited
        for (int i = 0; i < numberOfVertices; i++)
        {
            vertex[i]->setStatus(NOT_VISITED);
        }

        // Returns distance to destination vertex if it is reached
        return weight[dest];
    }
    // vector contains the distance from vertex[source] to all other vertices
    std::vector<int> distance(numberOfVertices);
    for (int i = 0; i < numberOfVertices; i++)
    {
        distance[i] = infinity;
    }
    distance[source] = 0;

    // Priority queue to store vertices and their distances
    Heap distancePriorityQueue;
    distancePriorityQueue.enqueue(new HNode(source, 0));
    while (!distancePriorityQueue.empty())
    {
        // Extract the vertex with the minimum distance from the priority queue
        HNode *current = distancePriorityQueue.dequeue();
        int vertexIndex = current->index;

        // If the vertex has been visited, skip it
        if (distance[vertexIndex] == infinity)
            continue;

        // Mark the vertex as visited
        distance[vertexIndex] = current->weight;
        vertex[vertexIndex]->setStatus(VISITED);
        delete current;

        // Explore neighbors
        for (int neighbor = 0; neighbor < numberOfVertices; ++neighbor)
        {
            // Relaxation step
            if (distance[neighbor] > (distance[vertexIndex] + distMatrix[vertexIndex][neighbor]))
            {
                distance[neighbor] = distance[vertexIndex] + distMatrix[vertexIndex][neighbor];
                distancePriorityQueue.enqueue(new HNode(neighbor, distance[neighbor]));
            }
        }
    }
    return distance[dest];
}
