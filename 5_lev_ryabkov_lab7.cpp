class Solution
{
public:
    int minimumDeviation(vector<int>& nums)
    {
        priority_queue<int> maxHeap; // O(N) по памяти — хранение всех элементов
        int minVal = INT_MAX; // O(1) по сложности, 4 байта

        for (int num : nums) // O(N log N) по сложности, O(N) по памяти (куча)
        {
            if (num % 2 == 1) num *= 2; // O(1) по сложности, 0 байт по памяти
            maxHeap.push(num); // O(log N) по сложности, O(N) по памяти (в сумме)
            minVal = min(minVal, num); // O(1) по сложности, 0 байт по памяти
        }

        int minDeviation = INT_MAX; // O(1) по сложности, 4 байта

        while (!maxHeap.empty() && maxHeap.top() % 2 == 0) // O(M log N) по сложности, 0 байт по памяти
        {
            int maxVal = maxHeap.top(); // O(1) по сложности, 4 байта
            maxHeap.pop(); // O(log N) по сложности, 0 байт по памяти

            minDeviation = min(minDeviation, maxVal - minVal); // O(1) по сложности, 0 байт по памяти
            maxVal /= 2; // O(1) по сложности, 0 байт по памяти
            maxHeap.push(maxVal); // O(log N) по сложности, O(N) по памяти (в сумме)
            minVal = min(minVal, maxVal); // O(1) по сложности, 0 байт по памяти
        }

        if (!maxHeap.empty()) // O(1) по сложности, 0 байт по памяти
        {
            minDeviation = min(minDeviation, maxHeap.top() - minVal); // O(1) по сложности, 0 байт по памяти
        }

        return minDeviation; // O(1) по сложности, 0 байт по памяти
    }
};

/*
Подсчёт памяти и сложности:
1. Куча (maxHeap):
   - O(N) по памяти для хранения всех элементов.
   - O(N log N) по сложности для вставки всех элементов.
2. Переменные minVal и minDeviation:
   - O(1) по сложности, 4 байта по памяти (на каждую переменную).
3. Основной цикл while:
   - O(M log N) по сложности, где M — количество операций деления максимального элемента.
   - O(1) по сложности для каждой операции в теле цикла.
   - 4 байта по памяти на локальную переменную maxVal.
4. Вспомогательные операции (умножение, деление, сравнение):
   - O(1) по сложности, 0 байт по памяти.

Итог:
- Сложность: O(N log N) в худшем случае.
- Память: O(N) для хранения элементов в куче + 8 байт на переменные.
*/