#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <deque>
#include <string>
#include <set>

// Simulate a data item created by processing a file
struct DataItem {
    int index;
    std::string content;

    DataItem(int idx, const std::string& cont) : index(idx), content(cont) {}
};

// Thread-safe queue for DataItem objects
class SafeQueue {
private:
    std::deque<DataItem> queue;
    std::mutex mtx;
    std::condition_variable cv;
    bool finished = false; // Flag to indicate no more items will be added

public:
    void push(const DataItem& item) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            queue.push_back(item);
        }
        cv.notify_one();
    }

    bool pop(DataItem& item) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&] { return !queue.empty() || finished; });
        if (!queue.empty()) {
            item = queue.front();
            queue.pop_front();
            return true;
        } else if (finished) {
            return false; // No more items will be added
        }
        return false;
    }

    void set_finished() {
        {
            std::lock_guard<std::mutex> lock(mtx);
            finished = true;
        }
        cv.notify_all();
    }
};

void producer_function(int index, SafeQueue& data_queue, std::mutex& output_mutex) {
    // Simulate file processing by creating a DataItem with some content
    std::string content = "Data from file " + std::to_string(index);
    DataItem item(index, content);

    // Simulate processing delay
    std::this_thread::sleep_for(std::chrono::milliseconds(100 * index));

    // Output message
    {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cout << "Producer " << index << " processed data: " << item.content << std::endl;
    }

    // Add the DataItem to the shared queue
    data_queue.push(item);
}

void consumer_function(SafeQueue& data_queue, std::mutex& output_mutex, int num_producers, std::vector<DataItem>& shared_data_items, std::mutex& data_items_mutex, std::set<std::pair<int, int>>& processed_pairs, std::mutex& pairs_mutex) {
    while (true) {
        DataItem item(0, "");
        if (!data_queue.pop(item)) {
            // No more items will be added and queue is empty
            break;
        }

        // Lock the shared data items vector
        {
            std::lock_guard<std::mutex> lock(data_items_mutex);

            // Compare the new item with all previous items
            for (const auto& prev_item : shared_data_items) {
                // To avoid duplicate processing, check if this pair has been processed
                int idx1 = prev_item.index;
                int idx2 = item.index;
                if (idx1 > idx2) std::swap(idx1, idx2);

                {
                    std::lock_guard<std::mutex> pairs_lock(pairs_mutex);
                    if (processed_pairs.find({idx1, idx2}) != processed_pairs.end()) {
                        continue; // Pair already processed
                    }
                    processed_pairs.insert({idx1, idx2});
                }

                // Simulate comparison
                std::string result = "Compared " + prev_item.content + " and " + item.content;

                // Output the result
                {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cout << "Consumer thread " << std::this_thread::get_id()
                              << " compared indices " << prev_item.index
                              << " and " << item.index << ": " << result << std::endl;
                }
            }

            // Add the new item to the list
            shared_data_items.push_back(item);
        }
    }
}

int main() {

    const int num_producers = 5;
    const int num_consumers = 3; // Number of consumer threads
    SafeQueue data_queue;
    std::mutex output_mutex;
    std::mutex data_items_mutex;
    std::vector<DataItem> shared_data_items;
    std::set<std::pair<int, int>> processed_pairs;
    std::mutex pairs_mutex;

    // Start producer threads
    std::vector<std::thread> producer_threads;
    for (int i = 1; i <= num_producers; ++i) {
        producer_threads.emplace_back(producer_function, i, std::ref(data_queue), std::ref(output_mutex));
    }

    // Start consumer threads
    std::vector<std::thread> consumer_threads;
    for (int i = 0; i < num_consumers; ++i) {
        consumer_threads.emplace_back(consumer_function, std::ref(data_queue), std::ref(output_mutex), num_producers, std::ref(shared_data_items), std::ref(data_items_mutex), std::ref(processed_pairs), std::ref(pairs_mutex));
    }

    // Wait for all producers to finish
    for (auto& t : producer_threads) {
        t.join();
    }

    // Signal to consumers that no more data will be added
    data_queue.set_finished();

    // Wait for consumers to finish
    for (auto& t : consumer_threads) {
        t.join();
    }

    return 0;
}
