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
        cv.wait(lock, [&] { return !queue.empty(); });
        item = queue.front();
        queue.pop_front();
        return true;
    }

    bool empty() {
        std::lock_guard<std::mutex> lock(mtx);
        return queue.empty();
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

void consumer_function(SafeQueue& data_queue, std::mutex& output_mutex, int num_producers) {
    std::vector<DataItem> data_items;
    int items_received = 0;

    while (items_received < num_producers) {
        DataItem item(0, "");
        data_queue.pop(item);

        // Compare the new item with all previous items
        for (const auto& prev_item : data_items) {
            // Simulate comparison
            std::string result = "Compared " + prev_item.content + " and " + item.content;

            // Output the result
            {
                std::lock_guard<std::mutex> lock(output_mutex);
                std::cout << "Consumer compared indices " << prev_item.index
                          << " and " << item.index << ": " << result << std::endl;
            }
        }

        // Add the new item to the list
        data_items.push_back(item);
        items_received++;
    }
}

int main() {
    const int num_producers = 5;
    SafeQueue data_queue;
    std::mutex output_mutex;

    // Start producer threads
    std::vector<std::thread> producer_threads;
    
    for (int i = 1; i <= num_producers; ++i) {
        producer_threads.emplace_back(producer_function, i, std::ref(data_queue), std::ref(output_mutex));
    }

    // Start consumer thread
    std::thread consumer_thread(consumer_function, std::ref(data_queue), std::ref(output_mutex), num_producers);

    // Wait for all producers to finish
    for (auto& t : producer_threads) {
        t.join();
    }

    // Wait for consumer to finish
    consumer_thread.join();

    return 0;
}
