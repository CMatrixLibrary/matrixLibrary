#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <functional>
#include <vector>
#include <future>

class ThreadPool {
private:
    std::vector<std::future<void>> tasks;

public:
    void addTask(std::function<void(void)> task) {
        tasks.push_back(std::async(std::launch::async, task));
    }

    void completeTasksAndStop() {
        for (auto& task : tasks) {
            task.wait();
        }
    }
};

#endif
