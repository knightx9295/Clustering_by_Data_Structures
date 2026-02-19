#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <list>
#include <algorithm>

using namespace std;

// struct for point
struct Point {
    double x, y;
};

// calculate Euclidean distance
double calculateDistance(const Point& p1, const Point& p2) {
    double dist = sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
    return round(dist * 10.0) / 10.0;  // Round to 1 decimal place
}

// Returns indices of closest points
vector<int> findClosestPointsIndices(const vector<vector<double>>& distanceMatrix, double& minValue) {
    int N = distanceMatrix.size();
    
    minValue = numeric_limits<double>::max();
    int minRow = -1, minCol = -1;
    
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            // Skip -1 values (merged groups)
            if (distanceMatrix[i][j] != -1 && distanceMatrix[i][j] < minValue) {
                minValue = distanceMatrix[i][j];
                minRow = i;
                minCol = j;
            }
        }
    }
    
    vector<int> indices;
    if (minRow != -1 && minCol != -1) {
        indices.push_back(minRow);
        indices.push_back(minCol);
    }
    
    return indices;
}

// Merge function as described in the algorithm
void mergeGroups(vector<vector<double>>& distanceMatrix, vector<list<int>>& groups, int idx1, int idx2) {
    int N = distanceMatrix.size();
    
    // Ensure idx1 < idx2 for consistent merging
    if (idx1 > idx2) {
        swap(idx1, idx2);
    }
    
    cout << "  Merging groups " << idx1+1 << " and " << idx2+1 << endl;
    
    // Merge the linked lists (O(1) operation)
    groups[idx1].splice(groups[idx1].end(), groups[idx2]);
    
    // Update distances for the new merged group
    for (int k = 0; k < N; k++) {
        if (k != idx1 && k != idx2) {
            // For upper triangular (i < j)
            if (idx1 < k) {
                // distanceMatrix[idx1][k] = min(distanceMatrix[idx1][k], distanceMatrix[idx2][k]);
                // But we need to be careful - idx2 might be < k or > k
                double dist1 = (idx2 < k) ? distanceMatrix[idx2][k] : distanceMatrix[k][idx2];
                double dist2 = distanceMatrix[idx1][k];
                
                if (dist1 != -1 && dist2 != -1) {
                    distanceMatrix[idx1][k] = min(dist1, dist2);
                } else if (dist1 != -1) {
                    distanceMatrix[idx1][k] = dist1;
                }
                // else keep dist2
                
                // Mark the old distance as invalid
                if (idx2 < k) {
                    distanceMatrix[idx2][k] = -1;
                } else {
                    distanceMatrix[k][idx2] = -1;
                }
            } 
            else if (k < idx1) {
                // distanceMatrix[k][idx1] = min(distanceMatrix[k][idx1], distanceMatrix[k][idx2]);
                double dist1 = distanceMatrix[k][idx1];
                double dist2 = distanceMatrix[k][idx2];
                
                if (dist1 != -1 && dist2 != -1) {
                    distanceMatrix[k][idx1] = min(dist1, dist2);
                } else if (dist2 != -1) {
                    distanceMatrix[k][idx1] = dist2;
                }
                // else keep dist1
                
                // Mark the old distance as invalid
                distanceMatrix[k][idx2] = -1;
            }
        }
    }
    
    // Mark the distance between the merged groups as -1
    distanceMatrix[idx1][idx2] = -1;
    
    // Mark the diagonal of merged group as -1
    distanceMatrix[idx2][idx2] = -1;
}

// Function to form groups until M groups remain
vector<list<int>> formGroups(vector<vector<double>>& distanceMatrix, const vector<Point>& points, int N, int M) {
    
    // Initialize groups: each point is its own group
    vector<list<int>> groups(N);
    for (int i = 0; i < N; i++) {
        groups[i].push_back(i);  // Store point index in the group
    }
    
    int currentGroups = N;
    int iteration = 1;
    
    cout << "\n=== Starting Group Formation ===\n";
    cout << "Initial: " << currentGroups << " groups\n";
    
    // Continue until we have M groups
    while (currentGroups > M) {
        cout << "\n--- Iteration " << iteration << " ---\n";
        cout << "Current groups: " << currentGroups << ", Target: " << M << endl;
        
        // Find the closest groups
        double minDist;
        vector<int> closestIndices = findClosestPointsIndices(distanceMatrix, minDist);
        
        if (closestIndices.empty()) {
            cout << "No valid groups to merge found.\n";
            break;
        }
        
        int idx1 = closestIndices[0];
        int idx2 = closestIndices[1];
        
        cout << "  Minimum distance: " << minDist 
             << " between groups " << idx1+1 << " and " << idx2+1 << endl;
        
        // Merge the groups
        mergeGroups(distanceMatrix, groups, idx1, idx2);
        
        // Decrease group count
        currentGroups--;
        
        // Print current state
        cout << "  Groups after merge: " << currentGroups << endl;
        
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
    
    // Fill the upper triangular part of the matrix
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dist = calculateDistance(points[i], points[j]);
            distanceMatrix[i][j] = dist;
        }
    }
    
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

int main() {
    string filename = "file.txt";  // Change to your filename
    vector<Point> points;
    int N, M;
    
    // Create the distance matrix
    vector<vector<double>> distanceMatrix = createDistanceMatrix(filename, points, N, M);
    
    cout << "========================================\n";
    cout << "        GROUPING ALGORITHM\n";
    cout << "========================================\n";
    cout << "Number of points (N): " << N << endl;
    cout << "Number of groups to form (M): " << M << endl;
    
    // Print initial state
    printPoints(points);
    cout << "\nInitial Distance Matrix:";
    printDistanceMatrix(distanceMatrix);
    
    // Form groups
    vector<list<int>> finalGroups = formGroups(distanceMatrix, points, N, M);
    
    // Print final groups
    printFinalGroups(finalGroups, points);
    
    // Print final distance matrix
    cout << "\nFinal Distance Matrix (after grouping):";
    printDistanceMatrix(distanceMatrix);
    
    return 0;
}