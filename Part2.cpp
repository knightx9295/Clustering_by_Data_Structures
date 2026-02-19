#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <list>
#include <algorithm>
#include <queue>

using namespace std;

// struct for point
struct Point {
    double x, y;
};

// Structure to store distance information in heap
struct HeapNode {
    double distance;
    int group1;
    int group2;
    
    // For min-heap comparison
    bool operator>(const HeapNode& other) const {
        return distance > other.distance;
    }
};

// MinHeap class implementation
class MinHeap {
private:
    vector<HeapNode> heap;
    
    void heapifyUp(int index) {
        while (index > 0) {
            int parent = (index - 1) / 2;
            if (heap[index].distance < heap[parent].distance) {
                swap(heap[index], heap[parent]);
                // Update location matrix
                int g1 = min(heap[index].group1, heap[index].group2);
                int g2 = max(heap[index].group1, heap[index].group2);
                int g1p = min(heap[parent].group1, heap[parent].group2);
                int g2p = max(heap[parent].group1, heap[parent].group2);
                swap(locationMatrix[g1][g2], locationMatrix[g1p][g2p]);
                index = parent;
            } else {
                break;
            }
        }
    }
    
    void heapifyDown(int index) {
        int size = heap.size();
        while (true) {
            int left = 2 * index + 1;
            int right = 2 * index + 2;
            int smallest = index;
            
            if (left < size && heap[left].distance < heap[smallest].distance) {
                smallest = left;
            }
            if (right < size && heap[right].distance < heap[smallest].distance) {
                smallest = right;
            }
            
            if (smallest != index) {
                swap(heap[index], heap[smallest]);
                // Update location matrix
                int g1 = min(heap[index].group1, heap[index].group2);
                int g2 = max(heap[index].group1, heap[index].group2);
                int g1s = min(heap[smallest].group1, heap[smallest].group2);
                int g2s = max(heap[smallest].group1, heap[smallest].group2);
                swap(locationMatrix[g1][g2], locationMatrix[g1s][g2s]);
                index = smallest;
            } else {
                break;
            }
        }
    }
    
public:
    // Constructor
    MinHeap(int N) {
        locationMatrix.resize(N, vector<int>(N, -1));
    }
    
    // Get size
    int size() const {
        return heap.size();
    }
    
    // Check if heap is empty
    bool empty() const {
        return heap.empty();
    }
    
    // Insert a new distance
    void insert(double distance, int group1, int group2) {
        // Ensure group1 < group2 for consistency
        if (group1 > group2) swap(group1, group2);
        
        HeapNode node;
        node.distance = distance;
        node.group1 = group1;
        node.group2 = group2;
        
        heap.push_back(node);
        int index = heap.size() - 1;
        locationMatrix[group1][group2] = index;
        
        heapifyUp(index);
    }
    
    // Extract minimum
    HeapNode extractMin() {
        if (heap.empty()) {
            return HeapNode{numeric_limits<double>::max(), -1, -1};
        }
        
        HeapNode minNode = heap[0];
        int lastIdx = heap.size() - 1;
        
        // Update location matrix
        int g1_min = min(minNode.group1, minNode.group2);
        int g2_min = max(minNode.group1, minNode.group2);
        locationMatrix[g1_min][g2_min] = -1;
        
        if (lastIdx > 0) {
            heap[0] = heap[lastIdx];
            int g1_last = min(heap[0].group1, heap[0].group2);
            int g2_last = max(heap[0].group1, heap[0].group2);
            locationMatrix[g1_last][g2_last] = 0;
        }
        
        heap.pop_back();
        
        if (!heap.empty()) {
            heapifyDown(0);
        }
        
        return minNode;
    }
    
    // Get heap index for a pair of groups
    int getHeapIndex(int group1, int group2) const {
        if (group1 > group2) swap(group1, group2);
        if (group1 >= locationMatrix.size() || group2 >= locationMatrix[group1].size()) {
            return -1;
        }
        return locationMatrix[group1][group2];
    }
    
    // Remove a specific distance
    void removeDistance(int group1, int group2) {
        if (group1 > group2) swap(group1, group2);
        
        int index = locationMatrix[group1][group2];
        if (index != -1 && index < heap.size()) {
            // Mark as very large value and heapify down
            heap[index].distance = numeric_limits<double>::max();
            heapifyDown(index);
            
            // If it's at the bottom, remove it
            if (index == heap.size() - 1) {
                heap.pop_back();
                locationMatrix[group1][group2] = -1;
            }
        }
    }
    
    // Update distance for a specific pair
    void updateDistance(int group1, int group2, double newDistance) {
        if (group1 > group2) swap(group1, group2);
        
        int index = locationMatrix[group1][group2];
        if (index != -1) {
            double oldDistance = heap[index].distance;
            heap[index].distance = newDistance;
            
            if (newDistance < oldDistance) {
                heapifyUp(index);
            } else {
                heapifyDown(index);
            }
        }
    }
    
    // Peek at minimum without extracting
    HeapNode peekMin() const {
        if (heap.empty()) {
            return HeapNode{numeric_limits<double>::max(), -1, -1};
        }
        return heap[0];
    }
    
    // Check if a distance exists in heap
    bool distanceExists(int group1, int group2) const {
        if (group1 > group2) swap(group1, group2);
        int index = getHeapIndex(group1, group2);
        return index != -1 && index < heap.size();
    }
    
    // Print heap
    void printHeap() const {
        cout << "Heap (size: " << heap.size() << "): ";
        for (size_t i = 0; i < min((size_t)10, heap.size()); i++) {
            cout << "[" << heap[i].group1+1 << "-" << heap[i].group2+1 
                 << ":" << heap[i].distance << "] ";
        }
        if (heap.size() > 10) cout << "...";
        cout << endl;
    }
    private:
    vector<vector<int>> locationMatrix;
    
};

// calculate Euclidean distance
double calculateDistance(const Point& p1, const Point& p2) {
    double dist = sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
    return round(dist * 10.0) / 10.0;  // Round to 1 decimal place
}

// Initialize min-heap with all pairwise distances
MinHeap initializeHeap(const vector<Point>& points, vector<vector<double>>& distanceMatrix) {
    int N = points.size();
    MinHeap heap(N);
    
    // Calculate all pairwise distances and insert into heap
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dist = calculateDistance(points[i], points[j]);
            distanceMatrix[i][j] = dist;
            heap.insert(dist, i, j);
        }
    }
    
    return heap;
}

// Find which group a point belongs to (for debugging)
int findGroup(const vector<list<int>>& groups, int pointIdx) {
    for (int i = 0; i < groups.size(); i++) {
        for (int p : groups[i]) {
            if (p == pointIdx) return i;
        }
    }
    return -1;
}

// Merge function using heap
void mergeGroupsWithHeap(vector<vector<double>>& distanceMatrix, vector<list<int>>& groups, MinHeap& heap,vector<int>& activeGroups,int idx1, int idx2) {
    int N = distanceMatrix.size();
    
    cout << "  Merging group " << idx1+1 << " (size: " << groups[idx1].size() 
         << ") and group " << idx2+1 << " (size: " << groups[idx2].size() << ")" << endl;
    
    // Merge the linked lists (O(1) operation)
    groups[idx1].splice(groups[idx1].end(), groups[idx2]);
    
    // Remove all distances involving group idx2 from heap
    for (int k = 0; k < N; k++) {
        if (k != idx2) {
            int g1 = min(idx2, k);
            int g2 = max(idx2, k);
            heap.removeDistance(g1, g2);
        }
    }
    
    // Remove distance between idx1 and idx2
    heap.removeDistance(min(idx1, idx2), max(idx1, idx2));
    
    // For each other active group, update distance to merged group
    for (int k = 0; k < N; k++) {
        if (k != idx1 && k != idx2 && !groups[k].empty()) {
            // Calculate minimum distance between merged group and group k
            double minDist = numeric_limits<double>::max();
            
            // Check all distances between points in merged group and points in group k
            for (int p1 : groups[idx1]) {
                for (int p2 : groups[k]) {
                    if (p1 < p2) {
                        minDist = min(minDist, distanceMatrix[p1][p2]);
                    } else if (p2 < p1) {
                        minDist = min(minDist, distanceMatrix[p2][p1]);
                    }
                }
            }
            
            if (minDist < numeric_limits<double>::max()) {
                // Update distance in matrix
                int g1 = min(idx1, k);
                int g2 = max(idx1, k);
                distanceMatrix[g1][g2] = minDist;
                
                // Update or insert in heap
                if (heap.distanceExists(g1, g2)) {
                    heap.updateDistance(g1, g2, minDist);
                } else {
                    heap.insert(minDist, g1, g2);
                }
            }
        }
    }
    
    // Mark group idx2 as inactive
    groups[idx2].clear();
    
    // Update active groups
    auto it = find(activeGroups.begin(), activeGroups.end(), idx2);
    if (it != activeGroups.end()) {
        activeGroups.erase(it);
    }
}

// Function to form groups until M groups remain using heap
vector<list<int>> formGroupsWithHeap(vector<vector<double>>& distanceMatrix, const vector<Point>& points, int N, int M) {
    
    // Initialize groups: each point is its own group
    vector<list<int>> groups(N);
    vector<int> activeGroups(N);
    for (int i = 0; i < N; i++) {
        groups[i].push_back(i);  // Store point index in the group
        activeGroups[i] = i;
    }
    
    // Initialize min-heap with all pairwise distances
    MinHeap heap = initializeHeap(points, distanceMatrix);
    
    int currentGroups = N;
    int iteration = 1;
    
    cout << "\n=== Starting Group Formation (Using Min-Heap) ===\n";
    cout << "Initial: " << currentGroups << " groups\n";
    cout << "Heap size: " << heap.size() << " distances\n";
    
    // Continue until we have M groups
    while (currentGroups > M && !heap.empty()) {
        cout << "\n--- Iteration " << iteration << " ---\n";
        cout << "Current groups: " << currentGroups << ", Target: " << M << endl;
        
        // Extract minimum distance from heap
        HeapNode minNode = heap.extractMin();
        
        // Skip if groups are not both active
        int idx1 = minNode.group1;
        int idx2 = minNode.group2;
        
        if (idx1 == -1 || idx2 == -1 || 
            groups[idx1].empty() || groups[idx2].empty()) {
            continue;
        }
        
        double minDist = minNode.distance;
        
        cout << "  Minimum distance: " << minDist 
             << " between groups " << idx1+1 << " and " << idx2+1 << endl;
        
        // Merge the groups
        mergeGroupsWithHeap(distanceMatrix, groups, heap, activeGroups, idx1, idx2);
        
        // Decrease group count
        currentGroups--;
        
        // Print current state
        cout << "  Groups after merge: " << currentGroups << endl;
        cout << "  Heap size after merge: " << heap.size() << endl;
        
        iteration++;
    }
    
    // Collect only non-empty groups
    vector<list<int>> finalGroups;
    for (int i = 0; i < N; i++) {
        if (!groups[i].empty()) {
            finalGroups.push_back(groups[i]);
        }
    }
    
    return finalGroups;
}

// Function to print the final groups with points
void printFinalGroups(const vector<list<int>>& finalGroups, const vector<Point>& points) {
    cout << "\n=== FINAL GROUPS ===\n";
    
    for (size_t i = 0; i < finalGroups.size(); i++) {
        cout << "Group " << i + 1 << ": ";
        
        bool first = true;
        for (int pointIdx : finalGroups[i]) {
            if (!first) {
                cout << ", ";
            }
            cout << "(" << points[pointIdx].x << ", " << points[pointIdx].y << ")";
            first = false;
        }
        cout << endl;
    }
    cout << "Total groups: " << finalGroups.size() << endl;
}

// Input function to read points and create distance matrix
vector<vector<double>> createDistanceMatrix(const string& filename, vector<Point>& points, int& N, int& M) {
    ifstream inputFile(filename);
    
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }
    
    // Read N and M from first line
    inputFile >> N >> M;
    
    // Read points
    points.resize(N);
    for (int i = 0; i < N; i++) {
        inputFile >> points[i].x >> points[i].y;
    }
    
    inputFile.close();
    
    // Create N×N distance matrix initialized with 0.0
    vector<vector<double>> distanceMatrix(N, vector<double>(N, 0.0));
    
    return distanceMatrix;
}

// Function to print the distance matrix
void printDistanceMatrix(const vector<vector<double>>& distanceMatrix) {
    int N = distanceMatrix.size();
    
    cout << "\nDistance Matrix:\n";
    cout << "    ";
    for (int j = 0; j < N; j++) {
        cout << setw(8) << "P" << j+1;
    }
    cout << endl;
    
    for (int i = 0; i < N; i++) {
        cout << "P" << i+1 << ":";
        for (int j = 0; j < N; j++) {
            if (j < i) {
                cout << setw(8) << " ";
            } else {
                if (distanceMatrix[i][j] == -1) {
                    cout << setw(8) << "-1";
                } else if (distanceMatrix[i][j] == 0) {
                    cout << setw(8) << "0.0";
                } else {
                    cout << fixed << setprecision(1) << setw(8) << distanceMatrix[i][j];
                }
            }
        }
        cout << endl;
    }
}

// Function to print points
void printPoints(const vector<Point>& points) {
    cout << "\nPoints:\n";
    for (size_t i = 0; i < points.size(); i++) {
        cout << "P" << i+1 << ": (" << points[i].x << ", " << points[i].y << ")\n";
    }
}

// Main function
int main() {
    string filename = "file.txt";  // Change to your filename
    vector<Point> points;
    int N, M;
    
    // Create empty distance matrix (will be filled by initializeHeap)
    vector<vector<double>> distanceMatrix = createDistanceMatrix(filename, points, N, M);
    
    cout << "========================================\n";
    cout << "  EFFICIENT GROUPING ALGORITHM (MIN-HEAP)\n";
    cout << "========================================\n";
    cout << "Number of points (N): " << N << endl;
    cout << "Number of groups to form (M): " << M << endl;
    
    // Print initial state
    printPoints(points);
    
    // Form groups using heap
    vector<list<int>> finalGroups = formGroupsWithHeap(distanceMatrix, points, N, M);
    
    // Print final groups
    printFinalGroups(finalGroups, points);
    
    // Print final distance matrix
    cout << "\nFinal Distance Matrix (after grouping):";
    printDistanceMatrix(distanceMatrix);
    
    return 0;
}