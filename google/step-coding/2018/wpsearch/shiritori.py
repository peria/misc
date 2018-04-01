#!/usr/bin/python3

import wp
import json


class Shiritori(object):

    def __init__(self):
        self.collection = wp.WikipediaCollection("data/wp.db")
        self.index = wp.Index("data/index.db", self.collection)
        self.category = None
        self.initial = None
        self.used = set()
        self.limit = 10

    def reset(self):
        self.category = None
        self.initial = None
        self.used = set()
        self.limit = 10

    def initialize(self, query):
        words = self.index.search_category(query)
        if words:
            self.category = query
            print('Category:{0}'.format(self.category))
            # Output hints
            print([word.reading for word in words])
            return self.reply('{0} のカテゴリでしりとりしましょう。何から初めますか。'.format(query))
        else:
            return self.reply('そのカテゴリは知りません。'.format(query))

    def ask_categories(self):
        categories = self.index.ask_categories()
        print(categories)
        self.reset()
        return self.reply('カテゴリーを変えましょう。オススメ カテゴリーは {0} です。'.format(categories[0]))

    def answer(self, query):
        word = self.index.get_word(self.category, query)
        if not word:
            return False, self.reply('{0}は{1}ではありません。私の勝ちです'.format(
                query, self.category))
        if self.initial and word.reading[0] != self.initial:
            return False, self.reply('{0}は {1} から始まりません。私の勝ちです'.format(
                query, self.initial))
        if query in self.used:
            return False, self.reply('既に {0}は言われています。私の勝ちです'.format(
                query))
        self.used.add(query)
        print(word.id, word.reading)

        # Accept user's answer
        last = word.reading[-1]
        print('last: ', last)
        candidates = self.index.search(self.category, last, self.limit)
        for cand in candidates:
            if cand.id in self.used:
                continue
            self.used.add(cand.id)
            self.initial = cand.reading[-1]
            print([w.id for w in self.index.search(self.category, self.initial, 10)])
            if cand.reading[-1] == 'ン':
              return False, self.reply('{0}。 ン で終わったので私の負けです'.format(cand.id))

            return True, self.reply(cand.id)

        return False, self.reply('もうありません。私の負けです')

    def thanks(self):
        if self.limit > 2:
            self.limit -= 1
        print('LIMIT {0}'.format(self.limit))
        return self.reply('どういたしまして')

    def reply(self, query='よく分かりません'):
        return json.dumps({
            'textToSpeech': query
        }, indent=2, separators=(',', ': '), ensure_ascii=False)

    @property
    def is_initialized(self):
        return bool(self.category)
