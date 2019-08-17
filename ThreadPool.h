#pragma once
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>

class ThreadPool {
private:
    std::mutex lock;
    std::condition_variable condVar;
    bool active = true;
    std::queue<std::function <void(void)>> tasks;
    std::vector<std::thread> threads;

public:
    ThreadPool(int numThread=std::thread::hardware_concurrency()) {
        if (numThread <= 0) numThread = 1;
        threads.reserve(numThread);
        for (int i = 0; i < numThread; ++i) {
            threads.emplace_back(std::bind(&ThreadPool::threadEntry, this));
        }
    }
    ~ThreadPool() {
        if (active) {
            completeTasksAndStop();
        }
    }

    void addTask(std::function<void(void)> task) {
        if (active) {
            std::unique_lock <std::mutex> l(lock);
            tasks.emplace(std::move(task));
            condVar.notify_one();
        }
    }

    void completeTasksAndStop() {
        {
            std::unique_lock <std::mutex> l(lock);
            active = false;
            condVar.notify_all();
        }
        for (auto& thread : threads) {
            thread.join();
        }
    }

private:
    void threadEntry() {
        std::function<void(void)> task;
        while (true) {
            {
                std::unique_lock<std::mutex> l(lock);

                while (active && tasks.empty()) condVar.wait(l);
                if (tasks.empty()) return;

                task = std::move(tasks.front());
                tasks.pop();
            }

            task();
        }
    }
};
